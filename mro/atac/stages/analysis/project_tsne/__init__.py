"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Obtain a TSNE projection of the reduced matrix
"""

import cellranger.analysis.pca as cr_pca
import cellranger.analysis.plsa as cr_plsa
import cellranger.analysis.lsa as cr_lsa
import cellranger.analysis.bhtsne as cr_tsne
import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io

from constants import ALLOWED_FACTORIZATIONS

import os
import tables
import numpy as np

__MRO__ = """
stage PROJECT_TSNE(
    in  h5       filtered_matrix,
    in  path     reduced_data,
    in  map      reduction_summary,
    in  int      tsne_perplexity,
    in  int      tsne_max_dims,
    in  int      tsne_input_pcs,
    in  float    tsne_theta,
    in  int      tsne_max_iter,
    in  int      tsne_stop_lying_iter,
    in  int      tsne_mom_switch_iter,
    in  int      random_seed,
    in  string[] factorization,
    out path     tsne,
    out map      tsne_summary,
    src py       "stages/analysis/project_tsne",
) split (
    in  string   method,
    in  int      tsne_dims,
) using (
    volatile = strict,
)
"""

def split(args):
    if args.filtered_matrix is None:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    if not os.path.exists(args.reduced_data):
        raise IOError('reduced data not found at {}'.format(args.reduced_data))

    if not set(args.factorization).issubset(ALLOWED_FACTORIZATIONS):
        raise ValueError('Invalid factorization provided')

    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix)

    reduction_summary = args.reduction_summary['h5']
    tsne_max_dims = analysis_constants.TSNE_N_COMPONENTS if args.tsne_max_dims is None else args.tsne_max_dims

    chunks = []
    for method in args.factorization:
        for tsne_dims in range(2, tsne_max_dims + 1):
            if method not in ALLOWED_FACTORIZATIONS:
                raise ValueError('invalid factorization method {}'.format(method))

            chunks.append({'method': method,
                           'tsne_dims': tsne_dims,
                           'transformed_matrix_h5': reduction_summary[method],
                           '__mem_gb': max(matrix_mem_gb, h5_constants.MIN_MEM_GB)})
    return {'chunks': chunks}

def main(args, outs):
    if args.filtered_matrix is None:
        return

    if not os.path.exists(outs.tsne):
        os.mkdir(outs.tsne)

    matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrix)

    if args.method == 'pca':
        transformed_matrix = cr_pca.load_pca_from_h5(args.transformed_matrix_h5).transformed_pca_matrix
    if args.method == 'lsa':
        lsa = cr_lsa.load_lsa_from_h5(args.transformed_matrix_h5)
        lsa = lsa._replace(transformed_lsa_matrix=lsa.transformed_lsa_matrix + 1e-120)
        # transform matrix to unit normed so that euclidean distance in new space is cosine distance in old space
        transformed_matrix = lsa.transformed_lsa_matrix / np.linalg.norm(lsa.transformed_lsa_matrix, axis=1, keepdims=True)
    if args.method == 'plsa':
        plsa = cr_plsa.load_plsa_from_h5(args.transformed_matrix_h5)
        plsa = plsa._replace(transformed_plsa_matrix=plsa.transformed_plsa_matrix + 1e-120)
        # transform matrix to unit normed so that euclidean distance in new space is cosine distance in old space
        transformed_matrix = plsa.transformed_plsa_matrix / np.linalg.norm(plsa.transformed_plsa_matrix, axis=1, keepdims=True)

    tsne_dims = args.tsne_dims
    tsne = cr_tsne.run_tsne(transformed_matrix,
                            key=str(tsne_dims),
                            tsne_dims=tsne_dims,
                            input_pcs=args.tsne_input_pcs,
                            perplexity=args.tsne_perplexity,
                            theta=args.tsne_theta,
                            max_iter=args.tsne_max_iter,
                            stop_lying_iter=args.tsne_stop_lying_iter,
                            mom_switch_iter=args.tsne_mom_switch_iter,
                            random_state=args.random_seed)

    filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
    _h5 = os.path.join(outs.tsne, args.method + '_tsne.h5')
    _csv = os.path.join(outs.tsne, args.method + '_tsne_csv')
    with tables.open_file(_h5, 'w', filters=filters) as f:
        cr_tsne.save_tsne_h5(tsne, f)

    cr_tsne.save_tsne_csv(tsne, matrix, _csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.filtered_matrix is None:
        outs.tsne = None
        outs.tsne_summary = {}
        return

    if not os.path.exists(outs.tsne):
        os.mkdir(outs.tsne)

    outs.tsne_summary = {'h5': {}, 'csv': {}}
    for method in args.factorization:
        # get all tsnes for a given method
        chunk_h5s = [os.path.join(chunk_out.tsne, method + '_tsne.h5')
                     for chunk_def, chunk_out in zip(chunk_defs, chunk_outs) if chunk_def.method == method]

        chunk_csv_dirs = [os.path.join(chunk_out.tsne, method + '_tsne_csv')
                          for chunk_def, chunk_out in zip(chunk_defs, chunk_outs) if chunk_def.method == method]

        analysis_io.combine_h5_files(chunk_h5s,
                                     os.path.join(outs.tsne, method + "_tsne.h5"),
                                     [analysis_constants.ANALYSIS_H5_TSNE_GROUP])

        for csv_dir in chunk_csv_dirs:
            cr_io.copytree(csv_dir,
                           os.path.join(outs.tsne, method + "_tsne_csv"),
                           allow_existing=True)

        outs.tsne_summary['h5'][method] = os.path.join(outs.tsne, method + "_tsne.h5")
        outs.tsne_summary['csv'][method] = os.path.join(outs.tsne, method + "_tsne_csv")
