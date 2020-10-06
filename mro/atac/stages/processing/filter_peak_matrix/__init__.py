"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Process and prune the peak matrix to relevant barcodes and peaks
"""

import os
import numpy as np
import martian
import cellranger.matrix as cr_matrix
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
import cellranger.atac.matrix as atac_matrix
import cellranger.library_constants as cr_lib_constants
import utils


__MRO__ = """
stage FILTER_PEAK_MATRIX(
    in  h5     raw_matrix,
    in  int    num_analysis_bcs,
    in  int    random_seed,
    in  csv    cell_barcodes,
    out h5     filtered_matrix,
    out path   filtered_matrix_mex,
    src py     "stages/processing/filter_peak_matrix",
)
"""

def split(args):
    if args.raw_matrix is None:
        return {'chunks': []}

    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.raw_matrix)

    res = {
        # multiplying factor to model proportion of observed barcodes that are cell barcodes
        '__mem_gb': max(round(np.ceil(1.3 * matrix_mem_gb)), h5_constants.MIN_MEM_GB),
    }
    return {'chunks': [], 'join': res}

def main(args, outs):
    martian.throw('No chunks defined.')

def prune(matrix, num_analysis_bcs=None, random_state=None):
    """Remove all cells that show no counts. If num_analysis_bcs is provided, it choses those number of barcodes.
    Finally, it returns a modified input matrix"""

    np.random.seed(0 if random_state is None else random_state)

    if matrix is None:
        return None

    if num_analysis_bcs:
        num_bcs = len(matrix.bcs)
        bc_indices = np.sort(np.random.choice(np.arange(num_bcs), size=min(num_analysis_bcs, num_bcs), replace=False))
        matrix = matrix.select_barcodes(bc_indices)

    nbcs = matrix.bcs_dim

    # keep barcode that have at least one peak left
    peaks_per_bc = (matrix.m > 0).sum(axis=0)
    keep_bcs = np.squeeze(np.array(peaks_per_bc != 0))

    matrix = matrix.select_barcodes(np.where(keep_bcs)[0])
    nbcs_new = matrix.bcs_dim
    martian.log_info("filtered out {} barcodes".format(nbcs - nbcs_new))

    return matrix

def join(args, outs, chunk_defs, chunk_outs):
    if args.raw_matrix is None:
        outs.filtered_matrix = None
        return

    # consume cell barcodes across all species and raise errors if not found
    if args.cell_barcodes is None:
        martian.exit("cell barcodes not provided")
    cell_barcodes = utils.load_cell_barcodes(args.cell_barcodes, with_species=True)

    # Read the peak matrix file and keep only cell barcodes
    # remove cell barcodes that were specified externally, such in reanalyzer,
    # which may not be present in raw matrix because they're missing from the fragments file
    present_cell_barcodes = {}
    peak_matrix = cr_matrix.CountMatrix.load_h5_file(args.raw_matrix)
    peak_matrix_bcs = set(peak_matrix.bcs)
    for species in cell_barcodes:
        present_cell_barcodes[species] = set()
        for bc in cell_barcodes[species]:
            if bc not in peak_matrix_bcs:
                martian.log_info("{} not found in the raw peak - bc matrix".format(bc))
            else:
                present_cell_barcodes[species].add(bc)

    peak_matrix = peak_matrix.filter_barcodes(present_cell_barcodes)
    if peak_matrix.features_dim == 0:
        martian.log_info("data has no peaks, skipping the clustering analysis")
        outs.filtered_matrix = None
        outs.filtered_matrix_mex = None
        return

    peak_matrix = prune(peak_matrix, num_analysis_bcs=args.num_analysis_bcs, random_state=args.random_seed)

    if peak_matrix.bcs_dim <= analysis_constants.MAX_N_CLUSTERS_DEFAULT:
        martian.log_info("Insufficient number of cell barcodes present after processing")
        outs.filtered_matrix = None
        outs.filtered_matrix_mex = None
        return

    if peak_matrix.features_dim < analysis_constants.MAX_N_CLUSTERS_DEFAULT:
        martian.log_info("Insufficient number of peaks present after processing")
        outs.filtered_matrix = None
        outs.filtered_matrix_mex = None
        return

    # save processed matrix
    peak_matrix.save_h5_file(outs.filtered_matrix, sw_version=martian.get_pipelines_version())
    if not os.path.exists(outs.filtered_matrix_mex):
        os.mkdir(outs.filtered_matrix_mex)
    atac_matrix.save_mex(peak_matrix, outs.filtered_matrix_mex,
                         cr_lib_constants.ATACSEQ_LIBRARY_TYPE,
                         sw_version=martian.get_pipelines_version())

    # save processed peaks
    # atac_feature_ref.save_features_bed_path(peak_matrix.feature_ref, outs.filtered_peaks)
