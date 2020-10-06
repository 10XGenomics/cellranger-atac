"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Cluster the cells using the reduced matrices from dimensionality reduction
"""

import numpy as np
import os

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.kmeans as cr_kmeans
import cellranger.analysis.pca as cr_pca
import cellranger.analysis.lsa as cr_lsa
import cellranger.analysis.plsa as cr_plsa
import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io

from constants import ALLOWED_FACTORIZATIONS, CLUSTER_FILE_HEAD

__MRO__ = """
stage CLUSTER_CELLS(
    in  path     reduced_data,
    in  map      reduction_summary,
    in  h5       filtered_matrix,
    in  string[] factorization,
    in  int      minclusters,
    in  int      maxclusters,
    in  int      num_dims,
    in  int      random_seed,
    out path     clustered_data,
    out map      clustering_summary,
    src py       "stages/analysis/cluster_cells",
)
"""

SAFETY_MEM_FACTOR = 1.4

def split(args):
    if args.filtered_matrix is None:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    if not os.path.exists(args.reduced_data):
        raise IOError('reduced data not found at {}'.format(args.reduced_data))

    if not set(args.factorization).issubset(ALLOWED_FACTORIZATIONS):
        raise ValueError('Invalid factorization provided')

    chunks = []
    min_clusters = args.minclusters if args.minclusters is not None else analysis_constants.MIN_N_CLUSTERS
    max_clusters = args.maxclusters if args.maxclusters is not None else analysis_constants.MAX_N_CLUSTERS_DEFAULT

    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    _, bcs_dim, _ = matrix_dims

    chunk_mem_gb = int(np.ceil(SAFETY_MEM_FACTOR * h5_constants.MIN_MEM_GB))
    reduction_summary = args.reduction_summary['h5']
    for n_clusters in xrange(min_clusters, max_clusters + 1):
        chunks.append({
            'transformed_matrix': reduction_summary,
            'n_clusters': n_clusters,
            '__mem_gb': chunk_mem_gb,
        })
    return {'chunks': chunks}

def main(args, outs):
    np.random.seed(0)

    if args.filtered_matrix is None:
        return

    if not os.path.exists(outs.clustered_data):
        cr_io.mkdir(outs.clustered_data)

    matrix_bcs = cr_matrix.CountMatrix.load_bcs_from_h5_file(args.filtered_matrix)
    for method in args.factorization:
        transformed_matrix = args.transformed_matrix[method]
        method_dir = os.path.join(outs.clustered_data, method)
        cr_io.mkdir(method_dir, allow_existing=True)
        file_head = CLUSTER_FILE_HEAD[method]
        _h5 = os.path.join(method_dir, file_head + ".h5")
        _csv = os.path.join(method_dir, file_head + "_csv")
        dr_mat = None

        if not os.path.exists(transformed_matrix):
            raise IOError('matrix does not exist')

        if method == 'pca':
            pca = cr_pca.load_pca_from_h5(transformed_matrix)
            dr_mat = pca.transformed_pca_matrix
        if method == 'lsa':
            lsa = cr_lsa.load_lsa_from_h5(transformed_matrix)
            lsa = lsa._replace(transformed_lsa_matrix=lsa.transformed_lsa_matrix + 1e-120)
            dr_mat = lsa.transformed_lsa_matrix / np.linalg.norm(lsa.transformed_lsa_matrix, axis=1, keepdims=True)
        if method == 'plsa':
            plsa = cr_plsa.load_plsa_from_h5(args.transformed_matrix[method])
            plsa = plsa._replace(transformed_plsa_matrix=plsa.transformed_plsa_matrix + 1e-120)
            dr_mat = plsa.transformed_plsa_matrix / np.linalg.norm(plsa.transformed_plsa_matrix, axis=1, keepdims=True)

        if args.num_dims is not None:
            if args.num_dims > dr_mat.shape[1]:
                raise ValueError('number of dimensions requested to use is larger than number of dimensions in data')
            dr_mat = dr_mat[:, np.arange(args.num_dims)]

        kmeans = cr_kmeans.run_kmeans(dr_mat, args.n_clusters, random_state=args.random_seed)
        with analysis_io.open_h5_for_writing(_h5) as f:
            cr_kmeans.save_kmeans_h5(f, args.n_clusters, kmeans)
        clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_KMEANS, args.n_clusters)
        cr_clustering.save_clustering_csv(_csv, clustering_key, kmeans.clusters, matrix_bcs)

def join(args, outs, chunk_defs, chunk_outs):
    if args.filtered_matrix is None:
        outs.clustered_data = None
        outs.clustering_summary = {}
        return

    if not os.path.exists(outs.clustered_data):
        cr_io.mkdir(outs.clustered_data)

    outs.clustering_summary = {'h5': {}, 'csv': {}}
    for method in args.factorization:
        chunk_h5s = [os.path.join(chunk_out.clustered_data, method, CLUSTER_FILE_HEAD[method] + ".h5") for chunk_out in chunk_outs]
        chunk_csv_dirs = [os.path.join(chunk_out.clustered_data, method, CLUSTER_FILE_HEAD[method] + "_csv") for chunk_out in chunk_outs]

        method_dir = os.path.join(outs.clustered_data, method)
        cr_io.mkdir(method_dir, allow_existing=True)
        analysis_io.combine_h5_files(chunk_h5s,
                                     os.path.join(method_dir, CLUSTER_FILE_HEAD[method] + ".h5"),
                                     [analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP,
                                      analysis_constants.ANALYSIS_H5_MAP_CLUSTERING[method]])

        for csv_dir in chunk_csv_dirs:
            cr_io.copytree(csv_dir,
                           os.path.join(method_dir, CLUSTER_FILE_HEAD[method] + "_csv"),
                           allow_existing=True)

        outs.clustering_summary['h5'][method] = os.path.join(method_dir, CLUSTER_FILE_HEAD[method] + ".h5")
        outs.clustering_summary['csv'][method] = os.path.join(method_dir, CLUSTER_FILE_HEAD[method] + "_csv")
