#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.analysis.stats as analysis_stats

import collections
from irlb import irlb
import numpy as np
import os
import tables

LSA = collections.namedtuple('LSA', ['transformed_lsa_matrix', 'components', 'variance_explained', 'dispersion', 'features_selected'])

def get_original_columns_used(cols_not_removed, cols_used_after_removal):
    """If a matrix  is subset down to only have columns indexed by cols_not_removed, and then is further subset to
    only contain cols_used_after removal, in that order, than this method returns the index of which columns in the old
    matrix correspond the the columns in the new matrix."""
    return [cols_not_removed[x] for x in cols_used_after_removal]

def run_lsa(matrix, lsa_features=None, lsa_bcs=None, n_lsa_components=None, random_state=None, discardPC=0, min_count_threshold=0):
    """ Run a LSA on the matrix using the IRLBA matrix factorization algorithm.  Prior to the LSA analysis, the counts
    are transformed by an inverse document frequency operation.

    If desired, only a subset of features (e.g. sample rows) can be selected for LSA analysis.  Each feature is ranked
    by its dispersion relative to other features that have a similar mean count.  The top `lsa_features` as ranked by
    this method will then be used for the LSA.

    One can also select to subset number of barcodes to use (e.g. sample columns), but in this case they are simply
    randomly sampled.

    Additionally one can choose to discard first N PCs (ranked by singular values/ variance explained). In this mode,
    the method automatically discovers N + n_lsa_components components

    Args:
        matrix (CountMatrix): The matrix to perform LSA on.
        lsa_features (int): Number of features to subset from matrix and use in LSA. The top lsa_features ranked by
                            dispersion are used
        lsa_bcs (int): Number of barcodes to randomly sample for the matrix.
        n_lsa_components (int): How many LSA components should be used.
        random_state (int): The seed for the RNG
        discardPC (int): number of components to discard
        min_count_threshold (int): The minimum sum of each row/column for that row/column to be passed to LSA
                                   (this filter is prior to any subsetting that occurs).
    Returns:
        A LSA object
    """
 
    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE
    np.random.seed(0)

    # Threshold the rows/columns of matrix, will throw error if an empty matrix results.
    thresholded_matrix, _, thresholded_features = matrix.select_axes_above_threshold(min_count_threshold)

    # If requested, we can subsample some of the barcodes to get a smaller matrix for LSA
    lsa_bc_indices = np.arange(thresholded_matrix.bcs_dim)
    if lsa_bcs is None:
        lsa_bcs = thresholded_matrix.bcs_dim
        lsa_bc_indices = np.arange(thresholded_matrix.bcs_dim)
    elif lsa_bcs < thresholded_matrix.bcs_dim:
        lsa_bc_indices = np.sort(np.random.choice(np.arange(thresholded_matrix.bcs_dim), size=lsa_bcs, replace=False))
    elif lsa_bcs > thresholded_matrix.bcs_dim:
        msg = ("You requested {} barcodes but the matrix after thresholding only "
               "included {}, so the smaller amount is being used.").format(lsa_bcs, thresholded_matrix.bcs_dim)
        print(msg)
        lsa_bcs = thresholded_matrix.bcs_dim
        lsa_bc_indices = np.arange(thresholded_matrix.bcs_dim)

    # If requested, select fewer features to use by selecting the features with highest normalized dispersion
    if lsa_features is None:
        lsa_features = thresholded_matrix.features_dim
    elif lsa_features > thresholded_matrix.features_dim:
        msg = ("You requested {} features but the matrix after thresholding only included {} features,"
               "so the smaller amount is being used.").format(lsa_features, thresholded_matrix.features_dim)
        print(msg)
        lsa_features = thresholded_matrix.features_dim
    # Calc mean and variance of counts after normalizing
    # But don't transform to log space, in order to preserve the mean-variance relationship
    m = analysis_stats.normalize_by_umi(thresholded_matrix)
    # Get mean and variance of rows
    (mu, var) = analysis_stats.summarize_columns(m.T)
    dispersion = analysis_stats.get_normalized_dispersion(mu.squeeze(), var.squeeze())  # TODO set number of bins?
    lsa_feature_indices = np.argsort(dispersion)[-lsa_features:]

    # Now determine how many components.
    if n_lsa_components is None:
        n_lsa_components = analysis_constants.LSA_N_COMPONENTS_DEFAULT

    # increment number of components if we discard PCs
    n_lsa_components += discardPC

    likely_matrix_rank = min(lsa_features, lsa_bcs)
    if likely_matrix_rank < n_lsa_components:
        print(("There are fewer nonzero features or barcodes ({}) than requested "
               "LSA components ({}); reducing the number of components.").format(likely_matrix_rank, n_lsa_components))
        n_lsa_components = likely_matrix_rank

    if (likely_matrix_rank * 0.5) <= float(n_lsa_components):
        print("Requested number of LSA components is large relative to the matrix size, an exact approach to matrix factorization may be faster.")

    # perform idf transform, which is suited for lsa
    lsa_mat = thresholded_matrix.select_barcodes(lsa_bc_indices).select_features(lsa_feature_indices)
    lsa_norm_mat = normalize_and_transpose(lsa_mat)
    (u, d, v, _, _) = irlb(lsa_norm_mat, n_lsa_components, random_state=random_state)

    # project the matrix to complete the transform: X --> X*v = u*d
    full_norm_mat = normalize_and_transpose(matrix)
    # Get a coordinate map so we know which columns in the old matrix correspond to columns in the new
    org_cols_used = get_original_columns_used(thresholded_features, lsa_feature_indices)
    transformed_irlba_matrix = full_norm_mat[:, org_cols_used].dot(v)[:, discardPC:]
    irlba_components = np.zeros((n_lsa_components - discardPC, matrix.features_dim))
    irlba_components[:, org_cols_used] = v.T[discardPC:, :]

    # calc proportion of variance explained
    variance_explained = np.square(d)[discardPC:] / np.sum(lsa_norm_mat.data**2)
    n_lsa_components = n_lsa_components - discardPC

    features_selected = np.array([f.id for f in matrix.feature_ref.feature_defs])[org_cols_used]

    # sanity check dimensions
    assert transformed_irlba_matrix.shape == (matrix.bcs_dim, n_lsa_components)
    assert irlba_components.shape == (n_lsa_components, matrix.features_dim)
    assert variance_explained.shape == (n_lsa_components,)

    return LSA(transformed_irlba_matrix, irlba_components, variance_explained, dispersion, features_selected)

def normalize_and_transpose(matrix):
    matrix.tocsc()

    m = analysis_stats.normalize_by_idf(matrix)

    # Transpose
    m = m.T

    return m

def get_irlb_mem_gb_from_matrix_dim(nonzero_entries):
    irlba_mem_gb = round(np.ceil(1.0 * nonzero_entries / analysis_constants.NUM_IRLB_MATRIX_ENTRIES_PER_MEM_GB)) + analysis_constants.IRLB_BASE_MEM_GB
    return max(h5_constants.MIN_MEM_GB, irlba_mem_gb)

def save_lsa_csv(lsa_map, matrix, base_dir):
    for n_components, lsa in lsa_map.iteritems():
        n_components_dir = os.path.join(base_dir, '%d_components' % n_components)
        cr_io.makedirs(n_components_dir, allow_existing=True)

        matrix_fn = os.path.join(n_components_dir, 'projection.csv')
        n_columns = lsa.transformed_lsa_matrix.shape[1]
        assert n_columns <= n_components
        matrix_header = ['Barcode'] + ['PC-%d' % (i + 1) for i in xrange(n_columns)]
        analysis_io.save_matrix_csv(matrix_fn, lsa.transformed_lsa_matrix, matrix_header,
                                    matrix.bcs)

        components_fn = os.path.join(n_components_dir, 'components.csv')
        components_header = ['PC'] + [f.id for f in matrix.feature_ref.feature_defs]
        analysis_io.save_matrix_csv(components_fn, lsa.components, components_header,
                                    range(1, n_components + 1))

        variance_fn = os.path.join(n_components_dir, 'variance.csv')
        variance_header = ['PC', 'Proportion.Variance.Explained']
        analysis_io.save_matrix_csv(variance_fn, lsa.variance_explained, variance_header,
                                    range(1, n_components + 1))

        dispersion_fn = os.path.join(n_components_dir, 'dispersion.csv')
        dispersion_header = ['Feature', 'Normalized.Dispersion']
        analysis_io.save_matrix_csv(dispersion_fn, lsa.dispersion, dispersion_header,
                                    [f.id for f in matrix.feature_ref.feature_defs])

        features_fn = os.path.join(n_components_dir, 'features_selected.csv')
        # TODO: there are two columns here, but only 1 entry in the header...BAD
        features_header = ['Feature']
        analysis_io.save_matrix_csv(features_fn, lsa.features_selected, features_header, range(1, len(lsa.features_selected) + 1))

def save_lsa_h5(lsa_map, f):
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_LSA_GROUP)
    for n_components, lsa in lsa_map.iteritems():
        analysis_io.save_h5(f, group, str(n_components), lsa)

def load_lsa_from_h5(filename):
    """ Load just the LSA info from an analysis h5 """
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_LSA_GROUP]
        # Just take the first LSA object, assuming we never have multiple
        for _, lsa in analysis_io.load_h5_iter(group, LSA):
            return lsa
