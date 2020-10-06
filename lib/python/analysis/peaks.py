"""
Methods for doing peakcalling.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
from __future__ import division, print_function, absolute_import

from collections import namedtuple
import numpy as np
from scipy import stats
from analysis.em_fitting import (
    MLE_dispersion, MLE_geometric_parameter, weighted_mean
)

# The default count threshold to use both for initializing the parameter estimation and
# to fall back on when there's not much count data.
DEFAULT_THRESHOLD = 15

MixedModelParams = namedtuple("MixedModelParams",
                              ["p_zero_noise",
                               "mean_bg_noise",
                               "alpha_bg_noise",
                               "p_signal",
                               "frac_zero_noise",
                               "frac_bg_noise"])


def simple_threshold_estimate(count_dict, cutoff=0.25):
    """Does a simple estimate for the threshold count:
    pick the threshold such that the lower noise region contains
    cutoff (default 25%) percentage of the total reads.
    """
    count_values = np.array(sorted([k for k in count_dict if k > DEFAULT_THRESHOLD]))
    count_weights = np.array([count_dict[k] for k in count_values])
    read_counts = count_values * count_weights
    cum_frac = read_counts.cumsum() / read_counts.sum()

    threshold = count_values[cum_frac <= cutoff][-1] \
        if any(cum_frac <= cutoff) else count_values[1]
    return threshold


def estimate_threshold(params, odds_ratio=20, maxval=1000):
    """Given a set of mixed model parameters and an odds ratio,
    find the lowest count that has at least that odds
    ratio of being from the signal rather than a background distribution.
    """
    xvals = np.arange(maxval)
    frac_signal = (1 - params.frac_zero_noise - params.frac_bg_noise)
    yvals_zero = stats.geom.pmf(xvals, params.p_zero_noise, loc=-1)
    yvals_bg = stats.nbinom.pmf(xvals, 1 / params.alpha_bg_noise,
                                1 / (1 + params.alpha_bg_noise * params.mean_bg_noise))
    yvals_sig = stats.geom.pmf(xvals, params.p_signal, loc=-1)
    yvals_zero = params.frac_zero_noise * yvals_zero
    yvals_bg = params.frac_bg_noise * yvals_bg
    yvals_sig = frac_signal * yvals_sig

    yvals_total_bg = yvals_zero + yvals_bg

    for val in xvals[::-1]:
        if yvals_sig[val] / yvals_total_bg[val] < odds_ratio:
            return val
    return None


def calculate_probability_distribution(count_dict, params):
    """Returns the normalized probability that each count comes from the signal or the
    noise distributions given some parameters for those.
    Counts where both PMFs are zero are returned with a (non-normalized)
    zero probability of coming from either."""
    count_values = np.array(sorted([k for k in count_dict]))

    zero_noise_pmf = stats.geom.pmf(count_values, params.p_zero_noise, loc=0)
    noise_pmf = stats.nbinom.pmf(count_values, 1 / params.alpha_bg_noise,
                                 1 / (1 + params.alpha_bg_noise * params.mean_bg_noise))
    signal_pmf = stats.geom.pmf(count_values, params.p_signal, loc=0)

    zero_noise_pmf *= params.frac_zero_noise
    noise_pmf *= params.frac_bg_noise
    signal_pmf *= (1 - params.frac_zero_noise - params.frac_bg_noise)

    output = np.vstack((zero_noise_pmf, noise_pmf, signal_pmf))
    mask = output.sum(axis=0) == 0
    output[:, mask] = 0.5
    output /= np.sum(output, axis=0)
    output[:, mask] = 0.0
    return output


def weighted_parameter_estimates(count_dict, prob_matrix):
    """Given N counts and a 3 by N probability matrix with the likelihood that
    each count is from signal or noise
    distributions, estimate the maximum likelihood parameters for the distributions
    """
    prob_zero_noise = prob_matrix[0, :]
    prob_bg_noise = prob_matrix[1, :]
    prob_signal = prob_matrix[2, :]

    zero_noise_mask = prob_zero_noise > 0
    bg_noise_mask = prob_bg_noise > 0
    signal_mask = prob_signal > 0

    count_values = np.array(sorted([k for k in count_dict]))
    count_weights = np.array([count_dict[k] for k in count_values])

    p_zero_noise = MLE_geometric_parameter(
        count_values[zero_noise_mask],
        count_weights[zero_noise_mask] *
        prob_zero_noise[zero_noise_mask])
    mean_bg_noise = weighted_mean(
        count_values[bg_noise_mask],
        count_weights[bg_noise_mask] *
        prob_bg_noise[bg_noise_mask])
    alpha_bg_noise = MLE_dispersion(
        count_values[bg_noise_mask],
        count_weights[bg_noise_mask] *
        prob_bg_noise[bg_noise_mask])
    p_signal = MLE_geometric_parameter(
        count_values[signal_mask],
        count_weights[signal_mask] *
        prob_signal[signal_mask])
    frac_zero_noise = (prob_zero_noise * count_weights).sum() / \
        (prob_matrix * count_weights).sum()
    frac_bg_noise = (prob_bg_noise * count_weights).sum() / \
        (prob_matrix * count_weights).sum()

    return MixedModelParams(p_zero_noise, mean_bg_noise, alpha_bg_noise, p_signal,
                            frac_zero_noise, frac_bg_noise)


def normalize_params(params):
    """Make sure that the larger of the two NB distributions is
    always annotated as signal"""
    if params.p_zero_noise < params.p_signal:
        return MixedModelParams(params.p_signal,
                                params.mean_bg_noise,
                                params.alpha_bg_noise,
                                params.p_zero_noise,
                                1 - (params.frac_zero_noise + params.frac_bg_noise),
                                params.frac_bg_noise)
    return params


def estimate_parameters(count_dict, epsilon=1e-6, maxiter=500, seed_params=None):
    """Uses an expectation-maximization method to estimate the mixed model parameters.
    An initial guess is made by segmenting the count data into three regions,
    then parameters are iteratively estimated from these."""
    threshold = simple_threshold_estimate(count_dict)

    count_values = np.array(sorted([k for k in count_dict]))

    prob_matrix = np.empty((3, len(count_values)), dtype=float)
    prob_matrix[0, :] = count_values <= DEFAULT_THRESHOLD
    prob_matrix[1, :] = (count_values <= threshold) & (count_values > DEFAULT_THRESHOLD)
    prob_matrix[2, :] = count_values > threshold

    params = weighted_parameter_estimates(count_dict, prob_matrix) if seed_params is None else seed_params
    last_params = None
    i = 0
    # pylint: disable=unsubscriptable-object
    while i < maxiter and (
            last_params is None or any(
        [abs(last_params[k] - params[k]) / params[k] > epsilon for k in
         range(len(params))])
    ):
        i += 1
        prob_matrix = calculate_probability_distribution(count_dict, params)
        last_params = params
        params = normalize_params(weighted_parameter_estimates(count_dict, prob_matrix))
    return params


def estimate_final_threshold(count_dict, odds_ratio=20):
    """Provided with a count dictionary, estimate the peak threshold for the count data.
    """
    if len(count_dict) < 30 or sum(count_dict.values()) < 1e5:
        return DEFAULT_THRESHOLD, None
    TRIES = 3
    tries = 0
    params = estimate_parameters(count_dict)

    # if the zero inflation mean is smaller than negative binom noise mean, swap and redo
    while 1 / params.p_zero_noise < params.mean_bg_noise and tries < TRIES:
        new_params = MixedModelParams(params.p_signal,
                                      1 / params.p_zero_noise,
                                      params.alpha_bg_noise,
                                      1 / params.mean_bg_noise,
                                      1 - (params.frac_zero_noise + params.frac_bg_noise),
                                      params.frac_bg_noise)
        params = estimate_parameters(count_dict, seed_params=new_params)
        tries += 1

    return estimate_threshold(params, odds_ratio), params
