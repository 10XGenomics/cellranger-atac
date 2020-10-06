"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Reduce the dimensions of the peak matrix
"""

import cellranger.analysis.pca as cr_pca
import cellranger.analysis.lsa as cr_lsa
import cellranger.analysis.plsa as cr_plsa
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
from constants import ALLOWED_FACTORIZATIONS
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import martian

import tables
import os

__MRO__ = """
stage REDUCE_DIMENSIONS(
    in  string[] factorization,
    in  int      num_dims,
    in  h5       filtered_matrix,
    in  int      num_bcs,
    in  int      num_features,
    in  int      random_seed,
    out path     reduced_data,
    out map      reduction_summary,
    src py       "stages/analysis/reduce_dimensions",
)
"""

def split(args):
    if args.filtered_matrix is None:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    if not os.path.exists(args.filtered_matrix):
        raise IOError('processed peak matrix not found at {}'.format(args.filtered_matrix))

    if not set(args.factorization).issubset(ALLOWED_FACTORIZATIONS):
        raise ValueError('Invalid factorization provided')

    # memory usage of just loading the full matrix
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    (features_dim, bcs_dim, nonzero_entries) = matrix_dims
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(bcs_dim, nonzero_entries)

    chunks = []
    for method in args.factorization:
        # mem alloc
        if method == "pca":
            method_mem_gb = cr_pca.get_irlb_mem_gb_from_matrix_dim(nonzero_entries)
        elif method == "lsa":
            method_mem_gb = cr_lsa.get_irlb_mem_gb_from_matrix_dim(nonzero_entries)
        elif method == "plsa":
            method_mem_gb = cr_plsa.get_plsa_mem_gb_from_matrix_dim(nonzero_entries)
        else:
            continue
        mem_gb = max(method_mem_gb + matrix_mem_gb, h5_constants.MIN_MEM_GB)

        # thread alloc
        if method == 'pca' or method == 'lsa':
            # HACK - give big jobs more threads in order to avoid overloading a node
            threads = cr_io.get_thread_request_from_mem_gb(mem_gb)
        if method == 'plsa':
            threads = -4

        chunks += [{
            'method': method,
            '__mem_gb': mem_gb,
            '__vmem_gb': mem_gb + 8,
            '__threads': threads,
        }]
    return {'chunks': chunks}

def main(args, outs):
    if args.filtered_matrix is None:
        return

    if not os.path.exists(outs.reduced_data):
        cr_io.mkdir(outs.reduced_data)

    method_dir = os.path.join(outs.reduced_data, args.method)
    cr_io.mkdir(method_dir, allow_existing=True)
    _h5 = os.path.join(method_dir, args.method + ".h5")
    _csv = os.path.join(method_dir, args.method + "_csv")

    matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrix)
    if args.method == "pca":
        pca = cr_pca.run_pca(matrix,
                             n_pca_components=args.num_dims,
                             pca_bcs=args.num_bcs,
                             pca_features=args.num_features,
                             random_state=args.random_seed)
        pca_key = args.num_dims if args.num_dims is not None else analysis_constants.PCA_N_COMPONENTS_DEFAULT
        pca_map = {pca_key: pca}
        filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
        with tables.open_file(_h5, 'w', filters=filters) as f:
            cr_pca.save_pca_h5(pca_map, f)
        cr_pca.save_pca_csv(pca_map, matrix, _csv)

    elif args.method == "lsa":
        lsa = cr_lsa.run_lsa(matrix,
                             n_lsa_components=args.num_dims,
                             lsa_bcs=args.num_bcs,
                             lsa_features=args.num_features,
                             random_state=args.random_seed)
        lsa_key = args.num_dims if args.num_dims is not None else analysis_constants.LSA_N_COMPONENTS_DEFAULT
        lsa_map = {lsa_key: lsa}
        filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
        with tables.open_file(_h5, 'w', filters=filters) as f:
            cr_lsa.save_lsa_h5(lsa_map, f)
        cr_lsa.save_lsa_csv(lsa_map, matrix, _csv)

    elif args.method == "plsa":
        # FIXME: temp directory from martian?
        plsa = cr_plsa.run_plsa(matrix,
                                outs.reduced_data,
                                n_plsa_components=args.num_dims,
                                threads=martian.get_threads_allocation(),
                                plsa_bcs=args.num_bcs,
                                plsa_features=args.num_features,
                                random_state=args.random_seed)
        plsa_key = args.num_dims if args.num_dims is not None else analysis_constants.PLSA_N_COMPONENTS_DEFAULT
        plsa_map = {plsa_key: plsa}
        filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
        with tables.open_file(_h5, 'w', filters=filters) as f:
            cr_plsa.save_plsa_h5(plsa_map, f)
        cr_plsa.save_plsa_csv(plsa_map, matrix, _csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.filtered_matrix is None:
        outs.reduced_data = None
        outs.reduction_summary = {}
        return

    if not os.path.exists(outs.reduced_data):
        cr_io.mkdir(outs.reduced_data)

    # copy chunk outs
    for chunk_out in chunk_outs:
        cr_io.copytree(chunk_out.reduced_data, outs.reduced_data, allow_existing=True)

    # Use final destinations to update summary
    outs.reduction_summary = {'h5': {}, 'csv': {}}
    for method in args.factorization:
        outs.reduction_summary['h5'][method] = os.path.join(outs.reduced_data, method, method + ".h5")
        outs.reduction_summary['csv'][method] = os.path.join(outs.reduced_data, method, method + "_csv")
