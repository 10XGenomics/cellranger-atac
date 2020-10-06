from __future__ import division
import numpy as np
import pandas as pd
import sys
from scipy.stats import t
from scipy.stats import fisher_exact
from scipy.sparse import csc_matrix
import numexpr as ne

import cellranger.analysis.diffexp as cr_diffexp

ALLOWED_NB2_MODELS = ['poisson', 'geometric', 'nb']

def sanitize(Y, X, alpha=None):
    """Ensure that the inputs to the numerical solver have correct dimensions and consistent with other inputs"""
    # sanity checks
    assert X.ndim == 2
    assert Y.ndim == 2

    if alpha is not None:
        assert alpha.ndim == 1

    # verify dims
    N = Y.shape[0]
    M = Y.shape[1]
    P = X.shape[1]
    assert X.shape[0] == N
    assert N >= P
    return N, M, P

def NB2_IRLS(Y, X, alpha, verbose=False):
    """
    Y: N x M matrix of cells x features (in scipys sparse form)
    alpha: 1 x M dispersion vector
    F: N X M offset numpy array which can be specified as a tuple of matrices (A, B)
       where F = A*B (this allows computing F on the fly for large datasets)
    X: design matrix with P covariates
    L: P x P prior numpy array that is positive semidefinite

    it solves Y[:, j] ~ NB(exp(X * B[:, j] + f[:, j]), 1*alpha[j], L[j]) for each feature j
    X is N X P, L is P X P and it returns B, a P x M matrix
    """
    # sanitize
    Y = csc_matrix(Y.astype('float'))
    N, M, P = sanitize(Y, X, alpha)

    # vars: O(MP + MP^2) ~ O(M)
    beta = np.zeros((P, M))
    varB = np.zeros((M, P, P))

    # criteria
    maxiter = 1500
    tol = 1e-6
    totaliters = 0

    # per GLM scratch: O(2PN + P^2 + P + M + 9N) ~ O(M) + O(N)
    Xt = X.transpose()
    Xp = np.zeros((P, N))
    XpX = np.zeros((P, P))
    Xpz = np.zeros(P)
    deviance = np.zeros(M)
    w = np.zeros(N)
    y = w.copy()
    wz = w.copy()
    mu = w.copy()
    eta = w.copy()
    y = w.copy()
    yr = w.copy()
    yp = w.copy()
    a0 = w.copy()

    # Run a GLM for each j
    for j in range(M):
        iters = 0
        delta = 1

        # init
        y = Y[:, j].A.squeeze()
        a0 = np.full(N, alpha[j])
        yp = y + 1 / a0
        yr = y + 1e-100

        # init
        mu = y.copy() + 0.1
        eta = np.log(mu)  # link function

        # common deviance terms that donot depend on estimates
        cd = 2 * np.sum(yr * np.log(yr) - yp * np.log(a0 * yp))

        # iterate
        while iters < maxiter and delta >= tol:
            iters += 1

            # working variable
            wz = 1. / (1. + a0 * mu)
            w = mu * wz  # mean-variance function
            wz = (eta - 1) * w + y * wz

            # solve the linear system
            Xp = Xt * w
            XpX = np.dot(Xp, X)
            Xpz = np.dot(Xt, wz)
            try:
                beta[:, j] = np.linalg.solve(XpX, Xpz)
            except:
                # occasionally the matrix can be singular, so we add a small regularizer
                beta[:, j] = np.linalg.solve(XpX + np.diag(np.full(XpX.shape[0], 1e-10)), Xpz)

            # update the means
            eta = np.dot(X, beta[:, j])  # linear
            mu = ne.evaluate('exp(eta)')  # inverse link

            # clamp for numerical stability
            mu[mu < 1e-300] = 1e-300
            mu[np.isinf(mu)] = 1e300

            # deviance that depends on the estimate
            old_dev = deviance[j]
            deviance[j] = cd + 2 * np.sum(yp * np.log(1 + a0 * mu) - y * eta, axis=0)
            delta = np.abs(deviance[j] - old_dev)

        # estimate varB
        varB[j, :, :] = np.linalg.inv(XpX)

        totaliters += iters

    if verbose:
        print "avg iters:", totaliters / M

    return beta, varB, deviance

def empirical_dispersion(y, threshold=1e-4):
    """Estimate empirical dispersion"""
    assert threshold > 0
    from cellranger.analysis.stats import summarize_columns
    (mu, var) = summarize_columns(y.T)
    alpha_est = np.maximum((var.squeeze() - mu.squeeze()) / (np.square(mu.squeeze() + 1e-100)), threshold)
    return alpha_est

def Wald_1x1_NB2(B, varB, df):
    """Perform Wald test for 1 vs rest grouping with NB2 regression"""
    var_a = varB[:, 0, 0].squeeze()
    var_b = varB[:, 1, 1].squeeze()
    cov_ab = varB[:, 0, 1].squeeze()

    # means and log2fc
    mean_a = np.nan_to_num(np.exp(B[0, :].squeeze()))
    mean_b = np.nan_to_num(np.exp(B[1, :].squeeze()))
    log2fc = (B[0, :].squeeze() - B[1, :].squeeze()) / np.log(2)

    # pval
    var_ab = np.sqrt((var_a + var_b - 2 * cov_ab))
    t_stat = (B[0, :].squeeze() - B[1, :].squeeze()) / var_ab
    p_values = t.sf(np.abs(t_stat), df)

    # adjust
    adj_p_values = p_values.copy()
    adj_p_values = cr_diffexp.adjust_pvalue_bh(p_values)
    return log2fc, t_stat, p_values, adj_p_values, mean_a, mean_b

def fisher_2x2(binarized_y, group_a, group_b, select_fisher, verbose=False):
    """Perform fisher's exact test on 2 x 2 contigency tables indicating presence or absence frequency of peak in cluster A vs B"""
    Na = len(group_a)
    Nb = len(group_b)

    # binarize for making the contingency tables
    na = binarized_y[:, group_a].sum(axis=1).A.squeeze()
    nc = Na - na
    nb = binarized_y[:, group_b].sum(axis=1).A.squeeze()
    nd = Nb - nb
    log2fc = np.zeros(len(select_fisher))
    p_values = log2fc.copy()

    if len(select_fisher) > 0 and verbose:
        print 'Running Fisher\'s exact test on {} features'.format(len(select_fisher))
        sys.stdout.flush()

    # run fisher test
    for i, nf in enumerate(select_fisher):
        _, p_values[i] = fisher_exact([[na[nf], nb[nf]], [nc[nf], nd[nf]]])
        # add the 1's and 2's to ensure no infs are produced. The log2fc for these cases aren't wholly meaningful anyway
        freq_ratio = (na[nf] + 1) * (Nb + 2) / (Na + 2) / (nb[nf] + 1)
        log2fc[i] = np.log2(freq_ratio)

    return log2fc, p_values


def NBGLM_differential_expression(y, group_a, group_b, model, test_params=None, verbose=False):
    """Perform NB GLM regression to determine differentially enriched features
    Args:
        y: scipy csc matrix
        group_a: indicator array of group 1
        group_b: indicator array of group 2
        model: selects the type of model to use. Available: Poisson, NB2, Geometric
        test_params: (optionally) encodes the covariates per cell and pre-determined dispersion
        verbose: controls verbosity of debug print out
    """
    p = 2
    q = 0
    N = y.shape[1]
    M = y.shape[0]
    Na = len(group_a)
    Nb = len(group_b)
    assert N >= p
    assert Na >= 1
    assert Nb >= 1
    assert N == (Na + Nb)

    if test_params is None:
        test_params = {}
    if 'cov' in test_params:
        X = np.zeros((N, p + len(test_params['cov'])))
        X[group_a, 0] = 1
        X[group_b, 1] = 1
        for cov in test_params['cov']:
            X[:, p] = cov.squeeze()
            p += 1
    else:
        X = np.zeros((N, p))
        X[group_a, 0] = 1
        X[group_b, 1] = 1

    # init alpha
    # NOTE: Can use shrinkage estimators on empirical MoMs
    if model == 'nb':
        alpha_est = test_params.get('alpha', empirical_dispersion(y))
        alpha_est[alpha_est <= 0] = 1e-6  # set zero dispersion to nonzero small value
    if model == 'poisson':
        alpha_est = np.full(M, 1e-4)
    if model == 'geometric':
        alpha_est = np.full(M, 1.)

    # estimate
    B, varB, alpha_est = NB2_IRLS(y.T, X, alpha_est, verbose=verbose)

    # get pvals:
    log2fc, t_stat, p_values, adj_p_values, mean_a, mean_b = Wald_1x1_NB2(B, varB, N - p - q)

    # adjust low counts by fisher's test
    binarized_y = y > 0
    select_fisher = np.where(adj_p_values > 0.49)[0]
    log2fc_fisher, p_values_fisher = fisher_2x2(binarized_y, group_a, group_b, select_fisher)
    log2fc[select_fisher] = log2fc_fisher
    p_values[select_fisher] = p_values_fisher

    # adjust pvals
    adj_p_values = p_values.copy()
    adj_p_values = cr_diffexp.adjust_pvalue_bh(p_values)

    de_result = pd.DataFrame({
        'tested': np.full(len(mean_a), True),
        'mean_a': mean_a,
        'mean_b': mean_b,
        'log2_fold_change': log2fc,
        'p_value': p_values,
        'adjusted_p_value': adj_p_values,
        't_stats': t_stat,
        'alpha': alpha_est,
    })

    return de_result


def run_differential_expression(mm, clusters, model, impute_rest=False, test_params=None, verbose=False):
    """ Compute differential expression for each cluster vs all other cells
        Args: matrix      - FeatureBCMatrix  :  feature expression data
              clusters    - np.array(int) :  1-based cluster labels"""
    if model not in ALLOWED_NB2_MODELS:
        raise ValueError("{} model is not supported".format(model))

    n_clusters = np.max(clusters)

    # Create a numpy array with 3*K columns;
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((mm.shape[0], 3 * n_clusters))

    de_result = None
    for cluster in xrange(1, 1 + n_clusters):

        print '*' * 80
        print 'Computing DE for cluster %d...' % cluster
        sys.stdout.flush()

        if n_clusters == 2 and cluster == 2 and impute_rest:
            # NOTE: under Wald test, the pvals should be symmetric
            # so one can derive the values for the "rest" group using previous de result
            all_de_results[:, 0 + 3 * (cluster - 1)] = de_result['mean_b']
            all_de_results[:, 1 + 3 * (cluster - 1)] = -1.0 * de_result['log2_fold_change']
            all_de_results[:, 2 + 3 * (cluster - 1)] = 1.0 - de_result['adjusted_p_value']
            return cr_diffexp.DIFFERENTIAL_EXPRESSION(all_de_results)

        in_cluster = clusters == cluster
        group_a = np.flatnonzero(in_cluster)
        group_b = np.flatnonzero(np.logical_not(in_cluster))

        # compute DE
        de_result = NBGLM_differential_expression(mm, group_a, group_b, model, test_params, verbose)
        all_de_results[:, 0 + 3 * (cluster - 1)] = de_result['mean_a']
        all_de_results[:, 1 + 3 * (cluster - 1)] = de_result['log2_fold_change']
        all_de_results[:, 2 + 3 * (cluster - 1)] = de_result['adjusted_p_value']

    return cr_diffexp.DIFFERENTIAL_EXPRESSION(all_de_results)
