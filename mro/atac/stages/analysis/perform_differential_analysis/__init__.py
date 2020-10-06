"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.diffexp as cr_diffexp
import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.matrix as cr_matrix
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
from constants import ALLOWED_FACTORIZATIONS
import analysis.diffexp_nb as nb2_diffexp
from utils import normalize_matrix, get_barcode_gc
import martian
from tools import ReferenceManager


__MRO__ = """
stage PERFORM_DIFFERENTIAL_ANALYSIS(
    in  bed      peaks,
    in  string   reference_path,
    in  h5       filtered_peak_bc_matrix,
    in  h5       filtered_tf_bc_matrix,
    in  string[] factorization,
    in  path     clustering,
    in  map      clustering_summary,
    out path     enrichment_analysis,
    out map      enrichment_analysis_summary,
    src py       "stages/analysis/perform_differential_analysis",
) split (
    in  string   method,
    in  string   clustering_key,
    in  int      cluster,
    out csv      tmp_diffexp,
)
"""


def split(args):
    ctg_mgr = ReferenceManager(args.reference_path)
    species = ctg_mgr.list_species()
    if args.filtered_peak_bc_matrix is None or len(species) > 1:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    chunks = []
    matrix_mem_gb = 0.
    if args.filtered_tf_bc_matrix is not None:
        matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_tf_bc_matrix) * 1.5
    matrix_mem_gb += cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_peak_bc_matrix)
    chunk_mem_gb = int(np.ceil(max(matrix_mem_gb, h5_constants.MIN_MEM_GB)))

    if not set(args.factorization).issubset(ALLOWED_FACTORIZATIONS):
        raise ValueError('Invalid factorization provided')

    # create a chunk for each method x clustering combo
    for method in args.factorization:
        clustering_h5 = args.clustering_summary['h5'][method]
        for key in SingleGenomeAnalysis.load_clustering_keys_from_h5(clustering_h5):
            clustering = SingleGenomeAnalysis.load_clustering_from_h5(clustering_h5, key)
            for cluster in set(clustering.clusters):
                chunks.append({
                    'method': method,
                    'clustering_key': key,
                    'cluster': cluster,
                    '__mem_gb': chunk_mem_gb,
                    '__vmem_gb': chunk_mem_gb + int(np.ceil(ctg_mgr.get_vmem_est())) + 1,
                    '__threads': 1,
                })

    return {'chunks': chunks, 'join': {'__mem_gb': 3}}

def main(args, outs):
    """Run this for each method x clustering key combination from split"""
    ctg_mgr = ReferenceManager(args.reference_path)
    species = ctg_mgr.list_species()
    if args.filtered_peak_bc_matrix is None or len(species) > 1:
        return

    # Load the peak-BC matrix and a clustering and perform DE
    peak_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_peak_bc_matrix)
    clustering_h5 = args.clustering_summary['h5'][args.method]
    clustering = SingleGenomeAnalysis.load_clustering_from_h5(clustering_h5, args.clustering_key)
    mask = clustering.clusters == args.cluster
    clustering.clusters[mask] = 1
    clustering.clusters[np.logical_not(mask)] = 2

    # find depth using peak matrix and normalize
    scale = np.array(peak_matrix.m.sum(axis=0)).squeeze()
    depth = (scale + 1) / np.median(scale)

    cov_peak = [np.log(depth)]
    diffexp_peak = nb2_diffexp.run_differential_expression(peak_matrix.m, clustering.clusters, model='poisson',
                                                           impute_rest=True, test_params={'cov': cov_peak}, verbose=True)

    # find empirical estimates of alpha
    tf_matrix = None
    diffexp_tf = None
    # do DE on tf-BC matrix
    if args.filtered_tf_bc_matrix is not None:
        tf_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_tf_bc_matrix)
        ntfmatrix = normalize_matrix(tf_matrix.m, scale)
        alpha_tf = nb2_diffexp.empirical_dispersion(ntfmatrix)
        barcode_GC = get_barcode_gc(args.reference_path, args.peaks, peak_matrix)
        cov_tf = [barcode_GC, np.log(depth)]
        diffexp_tf = nb2_diffexp.run_differential_expression(tf_matrix.m, clustering.clusters, model='nb', impute_rest=True,
                                                             test_params={'cov': cov_tf, 'alpha': alpha_tf}, verbose=True)

    # vstack
    diffexp = diffexp_peak if tf_matrix is None else cr_diffexp.DIFFERENTIAL_EXPRESSION(np.vstack([diffexp_peak.data, diffexp_tf.data]))

    # write out temp file
    np.savetxt(outs.tmp_diffexp, diffexp.data, delimiter=',')
    outs.enrichment_analysis = None
    outs.enrichment_analysis_summary = None

def join(args, outs, chunk_defs, chunk_outs):
    ctg_mgr = ReferenceManager(args.reference_path)
    species = ctg_mgr.list_species()
    if args.filtered_peak_bc_matrix is None or len(species) > 1:
        outs.enrichment_analysis = None
        outs.enrichment_analysis_summary = {}
        return

    peak_matrix_features = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(args.filtered_peak_bc_matrix)
    tf_matrix_features = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(args.filtered_tf_bc_matrix) if args.filtered_tf_bc_matrix is not None else None
    outs.enrichment_analysis_summary = {'h5': {}, 'csv': {}}
    # for each method, we merge h5 files and copy csv directories to one place
    cr_io.mkdir(outs.enrichment_analysis, allow_existing=True)
    for method in args.factorization:
        method_dir = os.path.join(outs.enrichment_analysis, method)
        cr_io.mkdir(method_dir, allow_existing=True)

        _h5 = os.path.join(method_dir, '{}_enrichment_h5.h5'.format(method))
        outs.enrichment_analysis_summary['h5'][method] = _h5
        chunk_h5s = []

        _csv = os.path.join(method_dir, '{}_enrichment_csv'.format(method))
        outs.enrichment_analysis_summary['csv'][method] = _csv
        diffexp_prefixes = [(fr.id, fr.name) for fr in peak_matrix_features.feature_defs]
        if args.filtered_tf_bc_matrix is not None:
            diffexp_prefixes += [(fr.id, fr.name) for fr in tf_matrix_features.feature_defs]

        clustering_h5 = args.clustering_summary['h5'][method]
        for key in SingleGenomeAnalysis.load_clustering_keys_from_h5(clustering_h5):

            chunk_outs_def_method_clustering = sorted([[chunk_out, chunk_def] for
                                                       chunk_out, chunk_def in zip(chunk_outs, chunk_defs)
                                                       if chunk_def.clustering_key == key], key=lambda x: x[1].cluster)
            chunk_outs_method_clustering = [c[0] for c in chunk_outs_def_method_clustering]

            # load 1 vs rest tests in sorted order of chunks and combine into one output per clustering
            diffexp = cr_diffexp.DIFFERENTIAL_EXPRESSION(np.hstack([np.loadtxt(com.tmp_diffexp, delimiter=',')[:, 0:3] for com in chunk_outs_method_clustering]))

            # write out h5
            chunk_h5 = martian.make_path('{}_enrichment_h5.h5'.format(key))
            with analysis_io.open_h5_for_writing(chunk_h5) as f:
                cr_diffexp.save_differential_expression_h5(f, key, diffexp)
            chunk_h5s += [chunk_h5]

            # write out csv
            cr_diffexp.save_differential_expression_csv_from_features(key, diffexp, diffexp_prefixes, _csv)

        analysis_io.combine_h5_files(chunk_h5s, _h5, [analysis_constants.ANALYSIS_H5_DIFFERENTIAL_EXPRESSION_GROUP,
                                                      analysis_constants.ANALYSIS_H5_MAP_DE[method]])
