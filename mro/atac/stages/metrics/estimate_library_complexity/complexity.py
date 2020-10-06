"""Function to compute library complexity"""

import numpy as np
from scipy import optimize

def interpolate_downsampled_counts(total_array, unique_array, ds_rate):
    '''Given a series of total counts and a series of unique counts, corresponding to downsampling rates between 0 and 1,
    find the downsampling interval for ds_rate and provided interpolated value for unique count'''

    assert ds_rate < 1.00001
    assert len(total_array) > 1
    assert len(total_array) == len(unique_array)

    # pylint: disable=undefined-loop-variable
    for num, t in enumerate(total_array[1:]):
        if ds_rate < t / total_array[-1]:
            break
    return (ds_rate * total_array[-1] - total_array[num]) / (total_array[num + 1] - total_array[num]) *\
           (unique_array[num + 1] - unique_array[num]) + unique_array[num]


def get_unique_and_total_fragments(dup_counts):
    """Sums up the unique and total fragments present from a dictionary of duplicate counts"""
    unique = 0
    total = 0
    for copy_num, read_count in dup_counts.iteritems():
        unique += read_count
        total += read_count * int(copy_num)
    return unique, total


def downsample_counts(dup_counts, sampling_rates=None):
    """Takes a dictionary of duplicate counts and a list of subsampling rates and produces numpy arrays of the
    unique and total fragments at each sampling rate, along with no downsampling
    """
    if sampling_rates is None:
        sampling_rates = [1.0]
    sampling_rates = [rate for rate in sampling_rates if rate <= 1.0]

    unique, total = get_unique_and_total_fragments(dup_counts)

    unique_array = np.empty(len(sampling_rates) + 1, dtype=float)
    total_array = np.empty(len(sampling_rates) + 1, dtype=float)
    unique_array[-1] = unique
    total_array[-1] = total

    def estimate_unique_loss(rate, copy_num):
        """Estimates the fraction of unique reads with a given number of copies that will become lost
        at a given downsampling rate
        """
        assert rate <= 1.0
        return 1 - np.power(1 - rate, copy_num)

    max_dup = int(max(dup_counts, key=int))
    duplicate_histogram = np.zeros(max_dup + 1)
    for copy_num, read_count in dup_counts.iteritems():
        duplicate_histogram[int(copy_num)] = read_count

    for i, rate in enumerate(sampling_rates):
        total_array[i] = total * rate
        unique_array[i] = (estimate_unique_loss(rate, np.arange(max_dup + 1)) * duplicate_histogram).sum()

    return unique_array, total_array


def mean_unique_from_total_dual(sampled, L1, L2, f):
    """Estimates the number of unique fragments found when sampling with replacement from a library with a,
    given number of unique fragments.  See http://lh3lh3.users.sourceforge.net/download/samtools.pdf section 1.1
    """
    U1 = L1 * (1 - np.exp(-(sampled * f) / L1))
    U2 = L2 * (1 - np.exp(-(sampled * (1 - f)) / L2))
    return U1 + U2


def mean_unique_from_total_safe(complexity_data):
    """Estimates the number of unique fragments found when sampling with replacement from a library with a
    given number of unique fragments.
    """
    if not np.isnan(complexity_data['fraction_pool1']):
        return mean_unique_from_total_dual(np.array(complexity_data['plot_total']),
                                           complexity_data['estimated_complexity1'],
                                           complexity_data['estimated_complexity2'],
                                           complexity_data['fraction_pool1'])
    else:
        return np.nan

def estimate_complexity_dual(sequenced, unique):
    """Uses least squares fitting to estimate the total number of unique fragments found in a library given
    sequenced versus unique fragment counts for overall and downsampled data.
    (Estimates relative populations of a 2-pool model)
    """
    def error(xdata, L1, L2, f):
        return mean_unique_from_total_dual(xdata, L1, L2, f)
    try:
        popt, cov = optimize.curve_fit(error, sequenced, unique, p0=[max(unique) / 2, max(unique) / 2, 0.0])
    except RuntimeError:
        return np.array([np.nan, np.nan, np.nan])

    # Handle fitting issues
    if popt[0] < 0 or popt[1] < 0 or popt[2] < 0 or popt[2] > 1:
        return np.array([np.nan, np.nan, np.nan])

    # if cv is too bad
    cv = np.sqrt(np.diag(cov)) / popt
    if cv[0] > 0.1 or cv[1] > 0.1:
        return np.array([np.nan, np.nan, np.nan])
    return popt


def estimate_library_complexity_safe(dup_data, method=None):
    """Wrapper around the 2-pool model"""
    sampling_rates = np.linspace(1e-5, 0.95, 100)

    if method == 'compressed':
        dup_hist = dup_data
    elif method == 'expanded':
        dup_hist = {int(k): sum(count for _, count in v) for k, v in dup_data.iteritems()}
    else:
        raise ValueError('Invalid format for duplicate data')

    failed_estimates = {'total': None,
                        'unique': None,
                        'estimated_complexity1': 0,
                        'estimated_complexity2': 0,
                        'fraction_pool1': np.nan,
                        'estimated_complexity': 0}

    if len(dup_hist) == 0:
        return failed_estimates
    unique, total = downsample_counts(dup_hist, sampling_rates)
    estimated_complexity1, estimated_complexity2, fraction_pool1 = estimate_complexity_dual(total, unique)

    # if the two-pool model failed, fall back on the single pool model
    if np.isnan(fraction_pool1):
        return {'total': total.tolist(),
                'unique': unique.tolist(),
                'estimated_complexity1': 0,
                'fraction_pool1': np.nan,
                'estimated_complexity2': 0,
                'estimated_complexity': 0}
    else:
        return {'total': total.tolist(),
                'unique': unique.tolist(),
                'estimated_complexity1': estimated_complexity1,
                'fraction_pool1': fraction_pool1,
                'estimated_complexity2': estimated_complexity2,
                'estimated_complexity': estimated_complexity1 + estimated_complexity2}
