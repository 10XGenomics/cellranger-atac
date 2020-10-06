#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
import martian
import math
import os

import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix
from constants import ALLOWED_FACTORIZATIONS, DEFAULT_FACTORIZATION

__MRO__ = """
stage ANALYZER_PREFLIGHT(
    in  bed      peaks,
    in  h5       filtered_peak_bc_matrix,
    in  string[] factorization,
    in  int      tsne_perplexity,
    in  int      random_seed,
    in  float    tsne_theta,
    in  int      tsne_mom_switch_iter,
    in  int      tsne_stop_lying_iter,
    in  int      tsne_max_dims,
    in  int      tsne_input_pcs,
    in  int      tsne_max_iter,
    in  int      max_clusters,
    in  int      num_components,
    in  int      num_dr_bcs,
    in  int      num_dr_features,
    in  float    neighbor_a,
    in  float    neighbor_b,
    in  int      graphclust_neighbors,
    src py       "stages/preflight/atac_analyzer",
)
"""

def option(arg, default):
    return arg if arg is not None else default

def main(args, outs):
    if args.filtered_peak_bc_matrix is None:
        martian.log_info("Filtered matrix not available for analysis")
        return

    if not os.path.exists(args.filtered_peak_bc_matrix):
        martian.exit("Filtered matrix {} does not exist".format(args.filtered_peak_bc_matrix))

    flt_matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_peak_bc_matrix)

    # check for empty matrix
    nonzero_entries = flt_matrix_dims[2]
    if nonzero_entries == 0:
        martian.exit("Peak-barcode matrix is empty - skipping analysis.")
        return

    total_peaks, total_bcs, _ = flt_matrix_dims

    # get parameters or their defaults
    factorization = option(args.factorization, [DEFAULT_FACTORIZATION])
    num_dr_bcs = option(args.num_dr_bcs, total_bcs)
    num_dr_features = option(args.num_dr_features, total_peaks)
    num_components = option(args.num_components, analysis_constants.LSA_N_COMPONENTS_DEFAULT)
    max_clusters = option(args.max_clusters, analysis_constants.MAX_N_CLUSTERS_DEFAULT)
    graphclust_neighbors = option(args.graphclust_neighbors, analysis_constants.GRAPHCLUST_NEIGHBORS_DEFAULT)
    neighbor_a = option(args.neighbor_a, analysis_constants.GRAPHCLUST_NEIGHBOR_A_DEFAULT)
    neighbor_b = option(args.neighbor_b, analysis_constants.GRAPHCLUST_NEIGHBOR_B_DEFAULT)
    tsne_pcs = option(args.tsne_input_pcs, num_components)
    tsne_max_iter = option(args.tsne_max_iter, analysis_constants.TSNE_MAX_ITER)
    tsne_mom_switch_iter = option(args.tsne_mom_switch_iter, analysis_constants.TSNE_MOM_SWITCH_ITER)
    tsne_stop_lying_iter = option(args.tsne_stop_lying_iter, analysis_constants.TSNE_STOP_LYING_ITER)
    tsne_theta = option(args.tsne_theta, analysis_constants.TSNE_THETA)
    tsne_max_dims = option(args.tsne_max_dims, analysis_constants.TSNE_N_COMPONENTS)
    tsne_perplexity = option(args.tsne_perplexity, analysis_constants.TSNE_DEFAULT_PERPLEXITY)

    # check constraints on parameters or their defaults
    # dim reduction
    if not set(factorization).issubset(ALLOWED_FACTORIZATIONS):
        martian.exit("factorization {} should be a subset of {}".format(','.join(factorization), ','.join(ALLOWED_FACTORIZATIONS)))
    if num_components > 20 and "plsa" in factorization:
        martian.exit("Recommend using fewer than 20 components while using PLSA for dimensionality reduction")
    if not (100 >= num_components >= 10):
        martian.exit("Number of components after dimensionality reduction must lie between 10 and 100")

    # dimensions
    if not (total_bcs >= num_dr_bcs >= max_clusters >= 2):
        martian.exit("""Parameters must satisfy total_bcs >= num_dr_bcs >= max_clusters >= 2. Possible causes:
        * the matrix has too few barcodes for analysis
        * you passed bad parameter values when calling 'cellranger-atac reanalyze'
        * you passed a barcode CSV file that's inconsistent with your matrix or analysis parameters""")

    if not (total_peaks >= num_dr_features >= num_components >= tsne_pcs >= 2):
        martian.exit("""Parameters must satisfy total_peaks >= num_dr_features >= num_components >= tsne_pcs >= 2. Possible causes:
        * the matrix has too few peaks for analysis (is your reference correct?)
        * you passed bad parameter values when calling 'cellranger-atac reanalyze'
        * you passed a peaks CSV file that's inconsistent with your matrix or analysis parameters""")

    # TSNE
    if not (tsne_max_iter >= tsne_mom_switch_iter >= 1 and tsne_max_iter >= tsne_stop_lying_iter >= 1):
        martian.exit("Parameters must satisfy tsne_max_iter >= tsne_mom_switch_iter >= 1 and tsne_max_iter >= tsne_stop_lying_iter >= 1.")
    if not (0 <= tsne_theta <= 1):
        martian.exit("Parameter tsne_theta must lie between 0 and 1.")
    if not (tsne_max_dims in [2, 3]):
        martian.exit("Parameter tsne_max_dims must be 2 or 3.")
    if not (1 <= tsne_perplexity <= 500):
        martian.exit("Parameter tsne_perplexity must lie between 1 and 500.")

    # clustering
    if not (max_clusters <= 50):
        martian.exit("Parameter max_clusters cannot be greater than 50.")
    if not (500 >= graphclust_neighbors >= 0):
        martian.exit("Parameter graphclust_neighbors must lie between 0 and 500.")
    if not (not math.isnan(neighbor_a) and not math.isinf(neighbor_a)):
        martian.exit("Parameter neighbor_a must be finite.")
    if not (neighbor_b >= 0):
        martian.exit("Parameter neighbor_b cannot be less than zero.")

    if not os.access(args.filtered_peak_bc_matrix, os.R_OK):
        martian.exit("Filtered matrix file is not readable, please check file permissions: %s" % args.filtered_peak_bc_matrix)
