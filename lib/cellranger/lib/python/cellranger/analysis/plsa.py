#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.analysis.stats as analysis_stats

import tenkit.log_subprocess as tk_subproc

import subprocess
import collections
import itertools
import numpy as np
import os
import tables

PLSA = collections.namedtuple('PLSA', ['transformed_plsa_matrix', 'components', 'variance_explained', 'dispersion', 'features_selected'])
PLSA_BINPATH = 'plsa'

def get_original_columns_used(cols_not_removed, cols_used_after_removal):
    """If a matrix  is subset down to only have columns indexed by cols_not_removed, and then is further subset to
    only contain cols_used_after removal, in that order, than this method returns the index of which columns in the old
    matrix correspond the the columns in the new matrix."""
    return [cols_not_removed[x] for x in cols_used_after_removal]

def run_plsa(matrix, temp_dir, plsa_features=None, plsa_bcs=None, n_plsa_components=None, random_state=None, threads=1, min_count_threshold=0):
    """ Run a PLSA on the matrix using the IRLBA matrix factorization algorithm.  Prior to the PLSA analysis, the
    matrix is not normalized at all.

    If desired, only a subset of features (e.g. sample rows) can be selected for PLSA analysis.  Each feature is ranked
    by its dispersion relative to other features that have a similar mean count.  The top `plsa_features` as ranked by
    this method will then be used for the PLSA.

    One *cannot* select to subset number of barcodes to use because of the intricacies of PLSA. It is still available as
    an optional input to match the API for lsa and pca subroutines included in this package.

    Args:
        matrix (CountMatrix): The matrix to perform PLSA on.
        plsa_features (int): Number of features to subset from matrix and use in PLSA. The top plsa_features ranked by
                            dispersion are used
        plsa_bcs (int): Number of barcodes to randomly sample for the matrix.
        n_plsa_components (int): How many PLSA components should be used.
        random_state (int): The seed for the RNG
        min_count_threshold (int): The minimum sum of each row/column for that row/column to be passed to PLSA
                                   (this filter is prior to any subsetting that occurs).
    Returns:
        A PLSA object
    """

    if not os.path.exists(temp_dir):
        raise Exception('Temporary directory does not exist. Need it to run plsa binary. Aborting..')

    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE
    np.random.seed(0)

    # Threshold the rows/columns of matrix, will throw error if an empty matrix results.
    thresholded_matrix, thresholded_bcs, thresholded_features = matrix.select_axes_above_threshold(min_count_threshold)

    # If requested, we can subsample some of the barcodes to get a smaller matrix for PLSA
    if plsa_bcs is not None:
        msg = "PLSA method does not allow subsetting barcodes"
        print(msg)
    plsa_bcs = thresholded_matrix.bcs_dim
    plsa_bc_indices = np.arange(thresholded_matrix.bcs_dim)

    # If requested, select fewer features to use by selecting the features with highest normalized dispersion
    if plsa_features is None:
        plsa_features = thresholded_matrix.features_dim
    elif plsa_features > thresholded_matrix.features_dim:
        msg = ("You requested {} features but the matrix after thresholding only included {} features,"
               "so the smaller amount is being used.").format(plsa_features, thresholded_matrix.features_dim)
        print(msg)
        plsa_features = thresholded_matrix.features_dim
    # Calc mean and variance of counts after normalizing
    # But don't transform to log space, in order to preserve the mean-variance relationship
    m = analysis_stats.normalize_by_umi(thresholded_matrix)
    # Get mean and variance of rows
    (mu, var) = analysis_stats.summarize_columns(m.T)
    dispersion = analysis_stats.get_normalized_dispersion(mu.squeeze(), var.squeeze())  # TODO set number of bins?
    plsa_feature_indices = np.argsort(dispersion)[-plsa_features:]

    # Now determine how many components.
    if n_plsa_components is None:
        n_plsa_components = analysis_constants.PLSA_N_COMPONENTS_DEFAULT

    likely_matrix_rank = min(plsa_features, plsa_bcs)
    if likely_matrix_rank < n_plsa_components:
        print(("There are fewer nonzero features or barcodes ({}) than requested "
               "PLSA components ({}); reducing the number of components.").format(likely_matrix_rank, n_plsa_components))
        n_plsa_components = likely_matrix_rank

    if (likely_matrix_rank * 0.5) <= float(n_plsa_components):
        print("Requested number of PLSA components is large relative to the matrix size, an exact approach to matrix factorization may be faster.")

    plsa_mat = thresholded_matrix.select_barcodes(plsa_bc_indices).select_features(plsa_feature_indices)

    # Write out sparse matrix without transforms
    # code picked up from save_mex
    plsa_mat.tocoo()
    out_matrix_fn = os.path.join(temp_dir, 'matrix.mtx')
    with open(out_matrix_fn, 'w') as stream:
        stream.write(np.compat.asbytes('%%MatrixMarket matrix {0} {1} {2}\n%%\n'.format('coordinate', 'integer', 'symmetry')))
        stream.write(np.compat.asbytes('%i %i %i\n' % (plsa_mat.m.shape[0], plsa_mat.m.shape[1], plsa_mat.m.nnz)))
        # write row, col, val in 1-based indexing
        for r, c, d in itertools.izip(plsa_mat.m.row + 1, plsa_mat.m.col + 1, plsa_mat.m.data):
            stream.write(np.compat.asbytes(("%i %i %i\n" % (r, c, d))))

    del plsa_mat

    # Run plsa module, reading in sparse matrix
    # Iters and tol are designed for 15PCs
    proc = tk_subproc.Popen([PLSA_BINPATH,
                             out_matrix_fn,
                             temp_dir,
                             '--topics', str(n_plsa_components),
                             '--iter', str(3000),
                             '--tol', str(0.002),
                             '--nt', str(threads)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_data, stderr_data = proc.communicate()
    if proc.returncode != 0:
        print stdout_data
        raise Exception("%s returned error code while running plsa binary %d: %s" % (proc, proc.returncode, stderr_data))

    # Read back data
    transformed_plsa_em_matrix_file = os.path.join(temp_dir, "transformed_matrix.csv")
    n_components_file = os.path.join(temp_dir, "components.csv")
    variance_explained_file = os.path.join(temp_dir, "topic_relevance.csv")
    org_rows_used = get_original_columns_used(thresholded_bcs, plsa_bc_indices)
    transformed_plsa_em_matrix = np.zeros((matrix.bcs_dim, n_plsa_components))
    transformed_plsa_em_matrix[org_rows_used, :] = np.genfromtxt(transformed_plsa_em_matrix_file, delimiter=",").astype('float64')
    org_cols_used = get_original_columns_used(thresholded_features, plsa_feature_indices)
    plsa_em_components = np.zeros((n_plsa_components, matrix.features_dim))
    plsa_em_components[:, org_cols_used] = np.genfromtxt(n_components_file, delimiter=",").astype('float64')
    variance_explained = np.genfromtxt(variance_explained_file, delimiter=",").astype('float64')

    # reorder components by variance explained as PLSA binary gives arbitrary order
    new_order = range(n_plsa_components)
    variance_explained, new_order = zip(*sorted(zip(variance_explained, new_order), reverse=True))
    variance_explained = np.array(variance_explained)
    plsa_em_components = plsa_em_components[new_order, :]
    transformed_plsa_em_matrix = transformed_plsa_em_matrix[:, new_order]

    # delete files
    cr_io.remove(transformed_plsa_em_matrix_file, allow_nonexisting=True)
    cr_io.remove(n_components_file, allow_nonexisting=True)
    cr_io.remove(variance_explained_file, allow_nonexisting=True)
    cr_io.remove(out_matrix_fn, allow_nonexisting=True)

    features_selected = np.array([f.id for f in matrix.feature_ref.feature_defs])[org_cols_used]

    # sanity check dimensions
    assert transformed_plsa_em_matrix.shape == (matrix.bcs_dim, n_plsa_components)
    assert plsa_em_components.shape == (n_plsa_components, matrix.features_dim)
    assert variance_explained.shape == (n_plsa_components,)

    return PLSA(transformed_plsa_em_matrix, plsa_em_components, variance_explained, dispersion, features_selected)

def get_plsa_mem_gb_from_matrix_dim(nonzero_entries):
    em_mem_gb = round(np.ceil(1.0 * nonzero_entries / analysis_constants.NUM_PLSA_EM_MATRIX_ENTRIES_PER_MEM_GB)) + analysis_constants.PLSA_EM_BASE_MEM_GB
    return max(h5_constants.MIN_MEM_GB, em_mem_gb)

def save_plsa_csv(plsa_map, matrix, base_dir):
    for n_components, plsa in plsa_map.iteritems():
        n_components_dir = os.path.join(base_dir, '%d_components' % n_components)
        cr_io.makedirs(n_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_components_dir, 'projection.csv')
        n_columns = plsa.transformed_plsa_matrix.shape[1]
        assert n_columns <= n_components
        matrix_header = ['Barcode'] + ['PC-%d' % (i + 1) for i in xrange(n_columns)]
        analysis_io.save_matrix_csv(matrix_fn, plsa.transformed_plsa_matrix, matrix_header,
                                    matrix.bcs)

        components_fn = os.path.join(n_components_dir, 'components.csv')
        components_header = ['PC'] + [f.id for f in matrix.feature_ref.feature_defs]
        analysis_io.save_matrix_csv(components_fn, plsa.components, components_header,
                                    range(1, n_components + 1))

        variance_fn = os.path.join(n_components_dir, 'variance.csv')
        variance_header = ['PC', 'Proportion.Variance.Explained']
        analysis_io.save_matrix_csv(variance_fn, plsa.variance_explained, variance_header,
                                    range(1, n_components + 1))

        dispersion_fn = os.path.join(n_components_dir, 'dispersion.csv')
        dispersion_header = ['Feature', 'Normalized.Dispersion']
        analysis_io.save_matrix_csv(dispersion_fn, plsa.dispersion, dispersion_header,
                                    [f.id for f in matrix.feature_ref.feature_defs])

        features_fn = os.path.join(n_components_dir, 'features_selected.csv')
        # TODO: there are two columns here, but only 1 entry in the header...BAD
        features_header = ['Feature']
        analysis_io.save_matrix_csv(features_fn, plsa.features_selected, features_header, range(1, len(plsa.features_selected) + 1))

def save_plsa_h5(plsa_map, f):
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_PLSA_GROUP)
    for n_components, plsa in plsa_map.iteritems():
        analysis_io.save_h5(f, group, str(n_components), plsa)

def load_plsa_from_h5(filename):
    """ Load just the PLSA info from an analysis h5 """
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_PLSA_GROUP]
        # Just take the first PLSA object, assuming we never have multiple
        for _, plsa in analysis_io.load_h5_iter(group, PLSA):
            return plsa
