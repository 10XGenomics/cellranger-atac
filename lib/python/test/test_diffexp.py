#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
# Test differential expression

import numpy as np
import scipy.sparse as sp_sparse
import tenkit.test as tk_test

import analysis.diffexp_nb as atac_de
from utils import normalize_matrix

def evaluate_de_res(de_res, de_genes, adj_p_cutoff=0.05):
    """ Compute sensitivty and PPV on DE calls in simulated data
        Args: de_res (pandas.DataFrame) - Result of sseq_differential_expression
              de_genes - Indices of genes that were truly DE (via simulation)
              adj_p_cutoff - Adjusted p-value cutoff
    """
    n_genes = de_res.shape[0]
    called_pos = np.logical_and(de_res['tested'], de_res['adjusted_p_value'] <= adj_p_cutoff)
    actual_pos = np.array([g in de_genes for g in xrange(n_genes)], dtype=bool)
    true_pos = np.logical_and(called_pos, actual_pos)
    sens = true_pos.sum().astype(float) / actual_pos.sum()
    ppv = true_pos.sum().astype(float) / called_pos.sum()

    return sens, ppv

def simulate_matrix():
    """ Simulate a gene-barcode matrix with 2 cell types but no biological variation.
    """
    # Number of cell types
    n_cell_types = 2
    # Total number of expressed genes
    n_genes = 10000
    # Log2 fold change of DE genes
    l2fc = 1.0
    # Number of DE genes per cell type
    n_de_genes = 100
    # Zipf / power law coefficient for simulated transcriptome
    txome_alpha = 0.7
    # Number of cells
    n_cells = 500
    # Mean/sd of transcripts-per-cell in log10 space
    n_transcripts_log10_mu = 5
    n_transcripts_log10_sd = 0.15
    # Conversion efficiency
    conv_eff = 0.1

    # Simulate the transcriptome via Zipf distribution
    base_txome = np.power(1 + np.arange(n_genes), - txome_alpha)
    base_txome = base_txome / np.sum(base_txome)

    # Generate the transcriptome for each cell-type by perturbing some genes
    de_genes = np.random.choice(n_genes, n_de_genes * n_cell_types, replace=False)
    txomes = np.zeros((n_genes, n_cell_types))
    for i in xrange(n_cell_types):
        txomes[:, i] = base_txome
        g = de_genes[(i * n_de_genes):(n_de_genes + i * n_de_genes)]
        txomes[g, i] = txomes[g, i] * np.power(2, l2fc)
        txomes[:, i] = txomes[:, i] / np.sum(txomes[:, i])

    # Assign cell types
    cell_type = np.random.choice(n_cell_types, n_cells, replace=True)

    # Generate the gene barcode matrix
    dense_mat = np.zeros((n_genes, n_cells), dtype=int)
    n = n_cells / n_cell_types

    for i in xrange(0, n_cells, n):
        # Simulate total number of transcripts for these cells
        n_transcripts = np.power(10, np.random.normal(n_transcripts_log10_mu, n_transcripts_log10_sd, n))

        # Simulate the transcriptomes for these cells
        type_counts = np.zeros((n, n_genes), dtype=int)
        for j in xrange(n):
            type_counts[j, :] = np.random.multinomial(n_transcripts[j], txomes[:, cell_type[i + j]])

        # Simulate mRNA capture
        obs_type_counts = np.random.binomial(type_counts, conv_eff)

        dense_mat[:, i:(i + n)] = np.transpose(obs_type_counts)

    mat = sp_sparse.csc_matrix(dense_mat)
    return mat, cell_type, de_genes


class TestDifferentialExpression(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_NB2_model(self):
        """Tests if the NB2 solver is correct by mixing together cell groups governed by a distinct N.B. distribution each
        but sharing the same dispersion parameter"""
        num_cells = [100000, 1000000]
        dispersion = 0.1
        prob1 = 0.2
        prob2 = 0.7

        for ncells in num_cells:
            split = int(ncells * 0.3)
            X = np.zeros((ncells, 2))

            # divide into 1/3 - 2/3 groups
            X[0:split, 0] = 1
            X[split:, 1] = 1

            Y = np.expand_dims(np.zeros(ncells), axis=1)
            Y[0:split, 0] = np.random.negative_binomial(int(1 / dispersion), prob1, size=split)
            Y[split:, 0] = np.random.negative_binomial(int(1 / dispersion), prob2, size=ncells - split)

            mu_expected = np.array([1 / dispersion * (1 - prob1) / prob1, 1 / dispersion * (1 - prob2) / prob2])

            BB, _, _ = atac_de.NB2_IRLS(Y, X, np.array([dispersion]), verbose=False)

            for expectation, estimation in zip(mu_expected, np.exp(BB)):
                self.assertApproxEqual(expectation, estimation, precision=1e-2)

    def test_simulated_gene_data(self):
        """ Test DE on a simulated gene expression matrix (w/ no biological variance) """
        np.random.seed(0)

        sim_mat, cell_type, sim_de = simulate_matrix()

        # get scale
        scale = np.array(sim_mat.sum(axis=0)).squeeze()
        depth = (scale + 1) / np.median(scale)
        cov = [np.log(depth)]

        # precompute distribution params
        ntfmatrix = normalize_matrix(sim_mat, scale)
        alpha = atac_de.empirical_dispersion(ntfmatrix)

        # sseq_params = cr_de.compute_sseq_params(sim_mat)
        # alpha = sseq_params['phi_g']

        de_res = atac_de.NBGLM_differential_expression(sim_mat, np.flatnonzero(cell_type == 0), np.flatnonzero(cell_type == 1),
                                                       model='nb', test_params={'cov': cov, 'alpha': alpha},
                                                       verbose=False)

        sensitivity, ppv = evaluate_de_res(de_res, sim_de)

        assert sensitivity >= 0.94
        assert ppv >= 0.94
