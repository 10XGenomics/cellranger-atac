#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
import martian
import numpy as np
import sys
import os
import h5py as h5

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.pca as cr_pca
import cellranger.analysis.lsa as cr_lsa
import cellranger.analysis.plsa as cr_plsa
import cellranger.matrix as cr_matrix
import cellranger.analysis.graphclust as cr_graphclust
import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
from cellranger.logperf import LogPerf
import cellranger.io as cr_io

from constants import ALLOWED_FACTORIZATIONS

__MRO__ = """
stage RUN_GRAPH_CLUSTERING(
    in  h5     matrix_h5           "Processed matrix",
    in  string[] factorization,
    in  path   reduced_data,
    in  int    num_neighbors       "Use this many neighbors",
    in  float  neighbor_a          "Use larger of (a+b*log10(n_cells) neighbors or num_neighbors",
    in  float  neighbor_b          "Use larger of (a+b*log10(n_cells) neighbors or num_neighbors",
    in  int    balltree_leaf_size,
    in  string similarity_type     "Type of similarity to use (nn or snn)",
    out h5     chunked_neighbors,
    out path   knn_clusters,
    out map    graph_clustering_summary,
    src py     "stages/analysis/run_graph_clustering",
) split using (
    in  string method,
    in  pickle neighbor_index,
    in  h5     submatrix,
    in  int    row_start,
    in  int    total_rows,
    in  int    k_nearest,
    in  h5     use_bcs,
)
"""

# 1e6 cells => ~64 chunks
NN_QUERIES_PER_CHUNK = 15000

DEFAULT_BALLTREE_LEAFSIZE = 40

# Memory usage in join, empirically determined
NN_ENTRIES_PER_MEM_GB = 10000000

# Unweighted nearest neighbor (boolean: is-nearest-neighbor)
NN_SIMILARITY = 'nn'

# Shared nearest neighbor (fraction of neighbors shared)
SNN_SIMILARITY = 'snn'

SIMILARITY_TYPES = [NN_SIMILARITY, SNN_SIMILARITY]

# TODO: Martian needs to provide a way to give split more memory.
# Workaround is mrp --overrides
def split(args):
    np.random.seed(0)

    if args.matrix_h5 is None:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    if not os.path.exists(args.reduced_data):
        raise IOError('reduced data not found at {}'.format(args.reduced_data))

    if not set(args.factorization).issubset(ALLOWED_FACTORIZATIONS):
        raise ValueError('Invalid factorization provided')

    if args.similarity_type not in SIMILARITY_TYPES:
        raise ValueError('Unsupported similarity type: %s. Must be one of: %s' % (args.similarity_type, ','.join(SIMILARITY_TYPES)))

    reduction_summary = args.reduction_summary['h5']

    method_dict = {}
    for method in args.factorization:
        method_dict[method] = {}

    with LogPerf('load'):
        for method in args.factorization:
            if method == 'pca':
                method_dict[method]['transformed_matrix'] = cr_pca.load_pca_from_h5(reduction_summary[method]).transformed_pca_matrix
            if method == 'lsa':
                method_dict[method]['transformed_matrix'] = cr_lsa.load_lsa_from_h5(reduction_summary[method]).transformed_lsa_matrix
            if method == 'plsa':
                method_dict[method]['transformed_matrix'] = cr_plsa.load_plsa_from_h5(reduction_summary[method]).transformed_plsa_matrix

    # Record indices of selected barcodes. All methods must use same barcodes
    use_bcs = np.arange(method_dict[args.factorization[0]]['transformed_matrix'].shape[0])
    use_bcs_path = martian.make_path('use_bcs.h5')
    cr_graphclust.save_ndarray_h5(use_bcs, use_bcs_path, 'use_bcs')

    # Build the nearest neighbor query index
    with LogPerf('nn_build'):
        for method in args.factorization:
            method_mat = method_dict[method]['transformed_matrix']
            # normalize for plsa/lsa so that standard euclidean distance in normalized space is cosine distance in original space
            if method in ['plsa', 'lsa']:
                method_mat = method_mat / np.linalg.norm(method_mat, axis=1, keepdims=True)
            balltree = cr_graphclust.build_neighbor_index(method_mat,
                                                          args.balltree_leaf_size or DEFAULT_BALLTREE_LEAFSIZE)
            method_dict[method]['neighbor_index'] = martian.make_path('neighbor_index_{}.pickle'.format(method))
            cr_graphclust.save_neighbor_index(balltree, method_dict[method]['neighbor_index'])

    # Compute the actual number of nearest neighbors we'll use
    given_num_neighbors = args.num_neighbors if args.num_neighbors is not None else analysis_constants.GRAPHCLUST_NEIGHBORS_DEFAULT
    given_neighbor_a = args.neighbor_a if args.neighbor_a is not None else analysis_constants.GRAPHCLUST_NEIGHBOR_A_DEFAULT
    given_neighbor_b = args.neighbor_b if args.neighbor_b is not None else analysis_constants.GRAPHCLUST_NEIGHBOR_B_DEFAULT

    # Take max of {num_neighbors, a + b*log10(n)}
    use_neighbors = int(max(given_num_neighbors, np.round(given_neighbor_a + given_neighbor_b * np.log10(len(use_bcs)))))

    # Clamp to [1, n - 1]
    num_neighbors = max(1, min(use_neighbors, len(use_bcs) - 1))
    print "Using %d neighbors" % num_neighbors

    # Divide the PCA matrix up into rows for NN queries
    with LogPerf('chunk_matrix'):
        chunks = []
        for method in args.factorization:
            method_mat = method_dict[method]['transformed_matrix']
            for row_start in xrange(0, method_mat.shape[0], NN_QUERIES_PER_CHUNK):
                row_end = min(row_start + NN_QUERIES_PER_CHUNK, method_mat.shape[0])

                # Write the submatrix to an h5 file
                submatrix_path = martian.make_path('{}_{}_submatrix.h5'.format(method, row_start))
                cr_graphclust.save_ndarray_h5(method_mat[row_start:row_end, :], submatrix_path, 'submatrix')

                chunks.append({
                    'method': method,
                    'neighbor_index': method_dict[method]['neighbor_index'],
                    'submatrix': submatrix_path,
                    'row_start': row_start,
                    'total_rows': method_mat.shape[0],
                    'k_nearest': num_neighbors,
                    'use_bcs': use_bcs_path,
                })

    if args.similarity_type == SNN_SIMILARITY:
        join_mem_gb = 64
        join_threads = 4  # Overallocate
    else:
        # Scale memory with size of nearest-neighbor adjacency matrix
        join_mem_gb = max(h5_constants.MIN_MEM_GB, int(np.ceil((num_neighbors * len(use_bcs)) / NN_ENTRIES_PER_MEM_GB)))
        # HACK: use more threads for bigger mem requests to avoid mem oversubscription on clusters that don't enforce it
        join_threads = cr_io.get_thread_request_from_mem_gb(join_mem_gb)

    return {
        'chunks': chunks,
        'join': {
            '__mem_gb': join_mem_gb,
            '__threads': join_threads,
        }}


def main(args, outs):
    np.random.seed(0)

    if args.matrix_h5 is None:
        outs.graph_clustering_summary = {}
        return

    with LogPerf('submatrix_load'):
        submatrix = cr_graphclust.load_ndarray_h5(args.submatrix, 'submatrix')

    with LogPerf('nn_idx_load'):
        balltree = cr_graphclust.load_neighbor_index(args.neighbor_index)

    with LogPerf('nn_query'):
        nn_matrix = cr_graphclust.compute_nearest_neighbors(submatrix, balltree, args.k_nearest, args.row_start)
        cr_graphclust.write_nearest_neighbors(nn_matrix, outs.chunked_neighbors)

def join(args, outs, chunk_defs, chunk_outs):
    if args.matrix_h5 is None:
        outs.graph_clustering_summary = {}
        return

    outs.graph_clustering_summary = {'h5': {}, 'csv': {}}
    # Merge the neighbor matrices
    for method in args.factorization:
        chunk_outs_def_method = [[chunk_out, chunk_def] for
                                 chunk_out, chunk_def in zip(chunk_outs, chunk_defs) if chunk_def.method == method]
        chunk_outs_method = [c[0] for c in chunk_outs_def_method]
        chunk_defs_method = [c[1] for c in chunk_outs_def_method]

        with LogPerf('merge_nn'):
            nn = cr_graphclust.merge_nearest_neighbors([chunk.chunked_neighbors for chunk in chunk_outs_method],
                                                       chunk_defs_method[0].total_rows)
        print 'nn\tnn_nodes\t%0.4f' % nn.shape[0]
        print 'nn\tnn_links\t%0.4f' % nn.nnz
        print 'nn\tnn_density\t%0.4f' % cr_graphclust.matrix_density(nn)
        sys.stdout.flush()

        matrix_bin = martian.make_path('matrix_{}.bin'.format(method))
        matrix_weights = martian.make_path('matrix_{}.weights'.format(method))
        louvain_out = martian.make_path('louvain_{}.out'.format(method))

        if args.similarity_type == 'snn':
            snn = cr_graphclust.compute_snn_matrix(nn, chunk_defs_method[0].k_nearest)

            print 'snn\tsnn_nodes\t%d' % snn.shape[0]
            print 'snn\tsnn_links\t%d' % (snn.nnz / 2)
            print 'snn\tsnn_density\t%0.4f' % ((snn.nnz) / float(snn.shape[0] * (snn.shape[0] - 1)))
            sys.stdout.flush()

            with LogPerf('convert'):
                cr_graphclust.pipe_weighted_edgelist_to_convert(snn, matrix_bin, matrix_weights)

            with LogPerf('louvain'):
                cr_graphclust.run_louvain_weighted_clustering(matrix_bin, matrix_weights, louvain_out)

        else:
            with LogPerf('tocoo'):
                nn = nn.tocoo(copy=False)

            with LogPerf('convert'):
                cr_graphclust.pipe_unweighted_edgelist_to_convert(nn, matrix_bin)

            with LogPerf('louvain'):
                cr_graphclust.run_louvain_unweighted_clustering(matrix_bin, louvain_out)

        with LogPerf('load_bcs'):
            barcodes = None
            with h5.File(args.matrix_h5, 'r') as f:
                group_name = f.keys()[0]
                barcodes = cr_matrix.CountMatrix.load_bcs_from_h5_group(f[group_name])

        use_bcs = cr_graphclust.load_ndarray_h5(chunk_defs_method[0].use_bcs, 'use_bcs')

        labels = cr_graphclust.load_louvain_results(len(barcodes), use_bcs, louvain_out)

        labels = cr_clustering.relabel_by_size(labels)

        # Save cluster results
        cr_io.mkdir(outs.knn_clusters, allow_existing=True)
        method_dir = os.path.join(outs.knn_clusters, method)
        cr_io.mkdir(method_dir, allow_existing=True)
        _h5 = os.path.join(method_dir, "clusters.h5")
        _csv = os.path.join(method_dir, "clusters_csv")
        with analysis_io.open_h5_for_writing(_h5) as f:
            cr_graphclust.save_graphclust_h5(f, labels)

        clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_GRAPHCLUST, 0)
        cr_clustering.save_clustering_csv(_csv, clustering_key, labels, barcodes)
        outs.graph_clustering_summary['h5'][method] = _h5
        outs.graph_clustering_summary['csv'][method] = _csv

    outs.chunked_neighbors = None
