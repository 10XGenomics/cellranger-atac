"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Tools for parameter estimation, used for EM model fitting.
"""

import numpy as np


def weighted_mean(counts, weights=None):
    """The mean of the input counts.  If weights are provided, returns the weighted mean."""
    if weights is None:
        return counts.mean()
    return (counts * weights).sum() / weights.sum()


def weighted_variance(counts, weights=None):
    """The variance of the input counts.  If weights are provided, returns the weighted variance."""
    if weights is None:
        return counts.var()
    sum_of_squares = (weighted_mean(counts, weights) - counts) ** 2
    return (sum_of_squares * weights).sum() / weights.sum()


def MOM_dispersion(counts, weights=None):
    """Simple --- and biased --- estimator of the dispersion parameter of a negative binomial distribution
    using the method of moments.
    """
    mu = weighted_mean(counts, weights)
    var = weighted_variance(counts, weights)
    alpha = (var - mu) / np.power(mu, 2)
    return max(1e-6, alpha)


def MLE_dispersion(counts, weights=None, maxiter=10000, epsilon=1e-6):
    """Borrowed from RiboPy (https://bitbucket.org/caethan/ribopy/src):
    Fixed point estimation of the dispersion (alpha) parameter for a negative binomial distribution.
    Provides the maximum likelihood estimate for alpha.
    """

    def G_mle(alpha, mu, counts, weights=None):
        """This is the fixed point for the maximum likelihood estimate of alpha (the NB dispersion parameter).

        alpha - float in (0, inf), the dispersion parameter
        mu - float, the fitted mean of the sample data
        counts - array of ints of length N, where N is the number of samples

        returns a single float that is a new estimate of alpha
        """
        v = np.arange(max(counts))
        subfunc_cumsum = np.cumsum((alpha * v) / (1 + alpha * v))

        indices = counts - 1
        indices[indices < 0] = 0

        return np.log(1 + alpha * mu) / (mu - weighted_mean(subfunc_cumsum[indices], weights))

    if (counts <= 0).all():
        return 1e-6
    alpha = MOM_dispersion(counts, weights)
    # If max(counts) is too big, this will be slow; fall back on method-of-moments
    if max(counts) > 1e8:
        return alpha

    mu = weighted_mean(counts, weights)
    var = weighted_variance(counts, weights)
    if var <= mu:
        return alpha
    for i in range(maxiter):
        new = G_mle(alpha, mu, counts, weights)
        if 2 * abs(new - alpha) / (new + alpha) < epsilon:
            break
        alpha = new
    return new


def MLE_geometric_parameter(counts, weights=None, exclude_zeros=True):
    """Returns the maximum likelihood estimate for a geometric fit to the data, with or without including zeros
    in the count data."""
    decrement = 1 if exclude_zeros else 0
    if weights is None:
        N = len(counts)
        k_total = sum(counts) - N * decrement
    else:
        N = sum(weights)
        k_total = sum(counts * weights) - N * decrement
    return N / (N + k_total)
