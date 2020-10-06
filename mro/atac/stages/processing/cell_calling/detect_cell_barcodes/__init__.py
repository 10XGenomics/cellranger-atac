"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Annotate barcodes which contain cells using total number of reads.

Uses an EM model to fit two negative binomial distributions to the reads-per-barcode data, one modeling the background
barcodes and one modeling the signal (cell-containing) barcodes.
"""
from __future__ import division, print_function

import json
from collections import namedtuple, Counter

import martian
import numpy as np
import pickle
from scipy import stats

from analysis.em_fitting import weighted_mean, MLE_dispersion
from barcodes import load_barcode_whitelist
from constants import MAXIMUM_CELLS_PER_SPECIES, NO_BARCODE, MINIMUM_COUNT, WHITELIST_CONTAM_RATE
from tenkit.stats import robust_divide
from tools import ReferenceManager, parsed_fragments_from_contig
from tools.regions import get_target_regions, fragment_overlaps_target

__MRO__ = """
stage DETECT_CELL_BARCODES(
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  string     barcode_whitelist,
    in  json       excluded_barcodes,
    in  map        force_cells,
    in  string     reference_path,
    in  bed        peaks,
    out csv        cell_barcodes,
    out csv        singlecell,
    out json       cell_calling_summary,
    src py         "stages/processing/cell_calling/detect_cell_barcodes",
) split (
    in  string     contig,
    out pickle     barcode_counts,
    out pickle     targeted_counts,
    out int        fragment_depth,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""
CELL_CALLING_THRESHOLD = 100000
MixedModelParams = namedtuple("MixedModelParams",
                              ["mu_noise", "alpha_noise", "mu_signal", "alpha_signal", "frac_noise"])


def split(args):
    if args.fragments is None:
        return {'chunks': [], 'join': {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    chunks = []
    for contig in all_contigs:
        chunks.append({'contig': contig,
                       '__mem_gb': 5})

    return {'chunks': chunks, 'join': {'__mem_gb': 5}}


def goodness_of_fit(count_data, params):
    """Goodness of fit metric using KS distance

    Count data is a simple array or list of counts
    """

    def mixed_nbinom_cdf(k):
        noise_cdf = stats.nbinom.cdf(k, 1 / params.alpha_noise, 1 / (1 + params.alpha_noise * params.mu_noise))
        signal_cdf = stats.nbinom.cdf(k, 1 / params.alpha_signal, 1 / (1 + params.alpha_signal * params.mu_signal))
        return params.frac_noise * noise_cdf + (1 - params.frac_noise) * signal_cdf

    gof, _ = stats.kstest(count_data, mixed_nbinom_cdf)
    return gof


def calculate_probability_distribution(count_dict, params):
    """Returns the normalized probability that each count comes from the signal or the noise distributions,
    given some parameters for those.  Counts where both PMFs are zero are returned with a (non-normalized)
    zero probability of coming from either."""
    count_values = np.array(sorted([k for k in count_dict]))

    noise_pmf = stats.nbinom.pmf(count_values,
                                 1 / params.alpha_noise,
                                 1 / (1 + params.alpha_noise * params.mu_noise))
    signal_pmf = stats.nbinom.pmf(count_values,
                                  1 / params.alpha_signal,
                                  1 / (1 + params.alpha_signal * params.mu_signal))

    output = np.vstack((params.frac_noise * noise_pmf, (1 - params.frac_noise) * signal_pmf))
    mask = output.sum(axis=0) == 0
    output[:, mask] = 0.5
    output /= np.sum(output, axis=0)
    output[:, mask] = 0.0
    return output


def weighted_parameter_estimates(count_dict, prob_matrix):
    """Given N counts and a 2 by N probability matrix with the likelihood that each count is from signal or noise
    distributions, estimate the maximum likelihood parameters for the distributions"""
    prob_noise = prob_matrix[0, :]
    prob_signal = prob_matrix[1, :]

    noise_mask = prob_noise > 0
    signal_mask = prob_signal > 0

    count_values = np.array(sorted([k for k in count_dict]))
    count_weights = np.array([count_dict[k] for k in count_values])

    mu_noise = weighted_mean(count_values[noise_mask], count_weights[noise_mask] * prob_noise[noise_mask])
    mu_signal = weighted_mean(count_values[signal_mask], count_weights[signal_mask] * prob_signal[signal_mask])
    alpha_noise = MLE_dispersion(count_values[noise_mask], count_weights[noise_mask] * prob_noise[noise_mask])
    alpha_signal = MLE_dispersion(count_values[signal_mask], count_weights[signal_mask] * prob_signal[signal_mask])
    frac_noise = (prob_noise * count_weights).sum() / (prob_matrix * count_weights).sum()
    return normalize_parameters(MixedModelParams(mu_noise, alpha_noise, mu_signal, alpha_signal, frac_noise))


def estimate_threshold(params, odds_ratio=20):
    """Estimate the separation threshold between signal and noise using the fits - look for the lowest count that
    gives an appropriate odds ratio in favor of signal
    """

    def ratio(k):
        signal = (1 - params.frac_noise) * stats.nbinom.pmf(k, 1 / params.alpha_signal,
                                                            1 / (1 + params.alpha_signal * params.mu_signal))
        noise = params.frac_noise * stats.nbinom.pmf(k, 1 / params.alpha_noise,
                                                     1 / (1 + params.alpha_noise * params.mu_noise))
        return signal / noise

    min_count = max(
        [
            int(params.mu_noise),
            stats.nbinom.ppf(0.001, 1 / params.alpha_signal, 1 / (1 + params.alpha_signal * params.mu_signal))
        ]
    )

    for test_count in np.arange(min_count, int(params.mu_signal) + 1):
        if ratio(test_count) >= odds_ratio:
            return test_count
    return test_count


def simple_threshold_estimate(count_dict, cutoff=0.01):
    """Does a simple estimate for the threshold count:  pick the threshold such that the lower noise region contains
    cutoff (default 1%) percentage of the total reads.
    """
    count_values = np.array(sorted([k for k in count_dict]))
    count_weights = np.array([count_dict[k] for k in count_values])
    read_counts = count_values * count_weights

    cum_frac = read_counts.cumsum() / read_counts.sum()

    threshold = count_values[cum_frac <= cutoff][-1] if any(cum_frac <= cutoff) else count_values[1]
    print(threshold, count_values[1], count_values[-2])
    threshold = max(threshold, count_values[1])
    threshold = min(threshold, count_values[-2])
    return threshold


def normalize_parameters(params):
    """Ensure that the larger of the two distributions is always the signal distribution."""
    if params.mu_noise > params.mu_signal:
        return MixedModelParams(params.mu_signal, params.alpha_signal,
                                params.mu_noise, params.alpha_noise,
                                1 - params.frac_noise)
    return params


def estimate_parameters(count_dict, epsilon=1e-6, maxiter=250, threshold=None):
    """Uses an expectation-maximization method to estimate the mixed model parameters."""
    if threshold is None:
        threshold = simple_threshold_estimate(count_dict)
    print('Initial threshold: {}'.format(threshold))

    count_values = np.array(sorted([k for k in count_dict]))

    prob_matrix = np.empty((2, len(count_values)), dtype=float)
    prob_matrix[0, :] = count_values <= threshold
    prob_matrix[1, :] = count_values > threshold
    params = weighted_parameter_estimates(count_dict, prob_matrix)
    print('Initial parameters')
    print(params)
    last_params = None
    i = 0
    while i < maxiter and (last_params is None or
                               any([abs(last_params[j] - params[j]) / params[j] > epsilon for j in
                                    range(len(params))])):
        i += 1
        prob_matrix = calculate_probability_distribution(count_dict, params)
        last_params = params
        params = weighted_parameter_estimates(count_dict, prob_matrix)
    print('Final parameters')
    print(params)
    return params


def compute_cell_index(species_list, cell_barcodes):
    """Produces a cell index dictionary that is keyed on barcode and links cell barcodes to per-species indices
    """
    cell_index = {}
    for species in species_list:
        bc_list = cell_barcodes.get(species, {}).keys()
        bc_list.sort()
        for index, barcode in enumerate(bc_list):
            label = "{species}_cell_{index}".format(**locals())
            if barcode not in cell_index:
                cell_index[barcode] = label
            else:
                cell_index[barcode] = "_".join([cell_index[barcode], label])
    return cell_index


def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()
    if args.fragments is None:
        outs.cell_barcodes = None
        outs.cell_calling_summary = None
        outs.singlecell = None
        return

    if args.excluded_barcodes is not None:
        with open(args.excluded_barcodes, 'r') as infile:
            excluded_barcodes = json.load(infile)
    else:
        excluded_barcodes = None

    # Merge the chunk inputs
    ref = ReferenceManager(args.reference_path)
    species_list = ref.list_species()

    barcode_counts_by_species = {species: Counter() for species in species_list}
    targeted_counts_by_species = {species: Counter() for species in species_list}
    fragment_depth = 0
    for chunk_in, chunk_out in zip(chunk_defs, chunk_outs):
        species = ref.species_from_contig(chunk_in.contig)
        with open(chunk_out.barcode_counts, 'r') as infile:
            barcode_counts_by_species[species] += pickle.load(infile)
        with open(chunk_out.targeted_counts, 'r') as infile:
            targeted_counts_by_species[species] += pickle.load(infile)
        fragment_depth += chunk_out.fragment_depth
    print('Total fragments across all chunks: {}'.format(fragment_depth))

    barcodes = list({bc
                     for species in species_list
                     for bc in barcode_counts_by_species[species]})
    non_excluded_barcodes = {species: [bc for bc in barcodes if bc not in excluded_barcodes[species]]
                             for species in species_list}
    print('Total barcodes observed: {}'.format(len(barcodes)))

    retained_counts = {}
    for species in species_list:
        if excluded_barcodes is None:
            retained_counts[species] = np.array([targeted_counts_by_species[species][bc] for bc in barcodes])
        else:
            retained_counts[species] = np.array([targeted_counts_by_species[species][bc] for bc in barcodes if
                                                 bc not in excluded_barcodes[species]])
            print('Barcodes excluded for species {}: {}'.format(species, len(excluded_barcodes[species])))
            print('Barcodes remaining for species {}: {}'.format(species, len(non_excluded_barcodes[species])))

    parameters = {}

    whitelist_length = len(load_barcode_whitelist(args.barcode_whitelist))
    count_shift = max(MINIMUM_COUNT, int(fragment_depth * WHITELIST_CONTAM_RATE / whitelist_length))
    print('Count shift for whitelist contamination: {}'.format(count_shift))

    for (species, count_data) in retained_counts.iteritems():
        print('Analyzing species {}'.format(species))
        # Subtract MINIMUM_COUNT from all counts to remove the effects of whitelist contamination
        shifted_data = count_data[count_data >= count_shift] - count_shift
        print('Number of barcodes analyzed: {}'.format(len(shifted_data)))
        count_dict = Counter(shifted_data)
        parameters[species] = {}

        forced_cell_count = None
        if args.force_cells is not None:
            if species in args.force_cells:
                forced_cell_count = int(args.force_cells[species])
            elif "default" in args.force_cells:
                forced_cell_count = int(args.force_cells["default"])
            if forced_cell_count > MAXIMUM_CELLS_PER_SPECIES:
                forced_cell_count = MAXIMUM_CELLS_PER_SPECIES
                martian.log_info(
                    'Attempted to force cells to {}.  Overriding to maximum allowed cells.'.format(forced_cell_count))

        # Initialize parameters to empty
        parameters[species]['noise_mean'] = None
        parameters[species]['noise_dispersion'] = None
        parameters[species]['signal_mean'] = None
        parameters[species]['signal_dispersion'] = None
        parameters[species]['fraction_noise'] = None
        parameters[species]['cell_threshold'] = None
        parameters[species]['goodness_of_fit'] = None
        parameters[species]['estimated_cells_present'] = 0

        # Corner case where FRIP is 0 because the number of peaks is tiny (fuzzer tests)
        if len(count_dict) < 10:
            parameters[species]['cells_detected'] = 0
            forced_cell_count = None
        elif forced_cell_count is None:
            print('Estimating parameters')
            fitted_params = estimate_parameters(count_dict)
            signal_threshold = estimate_threshold(fitted_params, CELL_CALLING_THRESHOLD) + count_shift
            print('Primary threshold: {}'.format(signal_threshold))
            parameters[species]['noise_mean'] = fitted_params.mu_noise
            parameters[species]['noise_dispersion'] = fitted_params.alpha_noise
            parameters[species]['signal_mean'] = fitted_params.mu_signal
            parameters[species]['signal_dispersion'] = fitted_params.alpha_signal
            parameters[species]['fraction_noise'] = fitted_params.frac_noise
            parameters[species]['cell_threshold'] = signal_threshold
            parameters[species]['goodness_of_fit'] = goodness_of_fit(shifted_data, fitted_params)
            called_cell_count = np.sum(count_data >= signal_threshold)
            parameters[species]['cells_detected'] = called_cell_count
            parameters[species]['estimated_cells_present'] = int((1 - fitted_params.frac_noise) * len(shifted_data))
            if called_cell_count > MAXIMUM_CELLS_PER_SPECIES:
                # Abort the model fitting and instead force cells to the maximum
                forced_cell_count = MAXIMUM_CELLS_PER_SPECIES

        if forced_cell_count is not None:
            print('Forcing cells to {}'.format(forced_cell_count))

            if forced_cell_count <= 0:
                raise ValueError("Force cells must be positive")
            else:
                adj_data = shifted_data[shifted_data > 0]
                print('Total barcodes considered for forcing cells: {}'.format(len(adj_data)))
                parameters[species]['cell_threshold'] = min(adj_data) if forced_cell_count >= len(adj_data) else \
                    sorted(adj_data, reverse=True)[forced_cell_count - 1]
                parameters[species]['cell_threshold'] += count_shift
                parameters[species]['cells_detected'] = np.sum(count_data >= parameters[species]['cell_threshold'])

    # For barnyard samples, mask out the noise distribution and re-fit to get cleaner separation
    if len(retained_counts) == 2 and (args.force_cells is None or not args.force_cells):
        print('Estimating secondary thresholds')
        sp1, sp2 = species_list

        sp1_threshold = -1 if parameters[sp1]['cell_threshold'] is not None else parameters[sp1]['cell_threshold']
        sp2_threshold = -1 if parameters[sp2]['cell_threshold'] is not None else parameters[sp2]['cell_threshold']

        if parameters[sp1]['cell_threshold'] is not None:
            sp1_counts = np.array([targeted_counts_by_species[sp1][bc] for bc in non_excluded_barcodes[sp1] if
                                   (targeted_counts_by_species[sp1][bc] > sp1_threshold) and
                                   (targeted_counts_by_species[sp2][bc] > sp2_threshold)])
            sp1_params = estimate_parameters(Counter(sp1_counts), threshold=sp1_threshold)
            if not np.isnan(sp1_params.frac_noise):
                parameters[sp1]['cell_threshold'] = max(sp1_threshold, estimate_threshold(sp1_params, 20))
            parameters[sp1]['cells_detected'] = np.sum(sp1_counts >= parameters[sp1]['cell_threshold'])
        else:
            parameters[sp1]['cells_detected'] = 0

        if parameters[sp2]['cell_threshold'] is not None:
            sp2_counts = np.array([targeted_counts_by_species[sp2][bc] for bc in non_excluded_barcodes[sp2] if
                                   (targeted_counts_by_species[sp1][bc] > sp1_threshold) and
                                   (targeted_counts_by_species[sp2][bc] > sp2_threshold)])
            sp2_params = estimate_parameters(Counter(sp2_counts), threshold=sp2_threshold)
            if not np.isnan(sp2_params.frac_noise):
                parameters[sp2]['cell_threshold'] = max(sp2_threshold, estimate_threshold(sp2_params, 20))
            parameters[sp2]['cells_detected'] = np.sum(sp2_counts >= parameters[sp2]['cell_threshold'])
        else:
            parameters[sp2]['cells_detected'] = 0

        print('Secondary threshold ({}): {}'.format(sp1, parameters[sp1]['cell_threshold']))
        print('Secondary threshold ({}): {}'.format(sp2, parameters[sp2]['cell_threshold']))

    print('Writing out cell barcodes')
    cell_barcodes = {}
    for (species, count_data) in retained_counts.iteritems():
        threshold = parameters[species]['cell_threshold']
        cell_barcodes[species] = {}
        print('Cell threshold for species {}: {}'.format(species, threshold))
        if threshold is not None:
            for count, barcode in zip(count_data, non_excluded_barcodes[species]):
                if count >= threshold:
                    print('{} - Total {}, Targeted {}, Count {}, Threshold {}'.format(
                        barcode, barcode_counts_by_species[species][barcode],
                        targeted_counts_by_species[species][barcode], count, threshold))
                    cell_barcodes[species][barcode] = count
        if len(cell_barcodes[species]) != parameters[species]['cells_detected']:
            print(len(cell_barcodes[species]), parameters[species]['cells_detected'])
            raise ValueError('Mismatch in called cells identified - failure in threshold setting')
        print('Selected {} barcodes of species {}'.format(len(cell_barcodes[species]), species))

    with open(outs.cell_barcodes, 'w') as outfile:
        # low mem reduce op to merge-sort bcs across species
        for species in cell_barcodes.keys():
            outfile.write(species + ",")
            outfile.write(",".join(cell_barcodes[species]) + "\n")

    cell_index = compute_cell_index(species_list, cell_barcodes)

    with open(outs.singlecell, 'w') as outfile:
        outfile.write("barcode,cell_id,")
        outfile.write(",".join(["is_{}_cell_barcode".format(species) for species in species_list]))
        if len(species_list) > 1:
            for species in species_list:
                outfile.write(",passed_filters_{}".format(species))
                outfile.write(",peak_region_fragments_{}".format(species))
        outfile.write("\n")
        for barcode in [NO_BARCODE] + sorted(barcodes):
            outfile.write("{},".format(barcode))
            outfile.write("{},".format(cell_index.get(barcode, "None")))
            values = [str(int(species in cell_barcodes and barcode in cell_barcodes[species]))
                      for species in species_list]
            outfile.write(",".join(values))
            if len(species_list) > 1:
                for species in species_list:
                    outfile.write(",{:d}".format(barcode_counts_by_species[species][barcode]))
                    outfile.write(",{:d}".format(targeted_counts_by_species[species][barcode]))
            outfile.write("\n")

    # process data into summary metrics
    summary_info = {}
    summary_info.update(generate_cell_calling_metrics(parameters, cell_barcodes))
    summary_info.update(generate_gb_metrics(cell_barcodes, excluded_barcodes))

    with open(outs.cell_calling_summary, 'w') as outfile:
        outfile.write(json.dumps(summary_info, indent=4))


def generate_gb_metrics(cell_barcodes, excluded_barcodes):
    """Estimate the GB doublet rate using cell count as a estimate of the total number
    of gel beads.
    """
    cell_barcode_set = set()
    for species, barcode_counts in cell_barcodes.iteritems():
        for barcode in barcode_counts:
            cell_barcode_set.add(barcode)
    total_cell_bcs = len(cell_barcode_set)

    excluded_gb_doublet_set = set()
    paired_gb_doublet_set = set()
    for species, exclusions in excluded_barcodes.iteritems():
        for excluded_bc, exclusion_data in exclusions.iteritems():
            reason, paired_bc = exclusion_data
            if reason != "gel_bead_doublet":
                continue
            assert excluded_bc not in cell_barcode_set
            excluded_gb_doublet_set.add(excluded_bc)
            paired_gb_doublet_set.add(paired_bc)

    estimated_doublet_gelbeads = len(excluded_gb_doublet_set)
    estimated_doublet_noncell_gelbeads = sum([bc not in cell_barcode_set for bc in paired_gb_doublet_set])

    estimated_total_gelbeads = total_cell_bcs + estimated_doublet_noncell_gelbeads
    doublet_rate = estimated_doublet_gelbeads / estimated_total_gelbeads

    metrics = {
        "estimated_gelbead_doublet_rate": doublet_rate,
        "fraction_gelbead_doublets_cells": 1 - (estimated_doublet_noncell_gelbeads / estimated_doublet_gelbeads),
    }
    return metrics


def generate_cell_calling_metrics(parameters, cell_barcodes):
    summary_info = {}
    species_list = parameters.keys()
    for species in species_list:
        key_suffix = "" if len(species_list) == 1 else "_{}".format(species)

        # Cell calling metrics
        summary_info["fitted_mean_noise{}".format(key_suffix)] = parameters[species]["noise_mean"]
        summary_info["fitted_dispersion_noise{}".format(key_suffix)] = parameters[species]["noise_dispersion"]
        summary_info["fitted_mean_signal{}".format(key_suffix)] = parameters[species]["signal_mean"]
        summary_info["fitted_dispersion_signal{}".format(key_suffix)] = parameters[species]["signal_dispersion"]
        summary_info["fraction_cell_calling_noise{}".format(key_suffix)] = parameters[species]["fraction_noise"]    
        summary_info["cell_threshold{}".format(key_suffix)] = parameters[species]["cell_threshold"]
        summary_info["goodness_of_fit{}".format(key_suffix)] = parameters[species]["goodness_of_fit"]
        summary_info["estimated_cells_present{}".format(key_suffix)] = parameters[species]["estimated_cells_present"]

        summary_info["annotated_cells{}".format(key_suffix)] = len(cell_barcodes[species])
        summary_info["estimated_fraction_cells_annotated{}".format(key_suffix)] = \
            robust_divide(len(cell_barcodes[species]), parameters[species]["estimated_cells_present"])

    summary_info["cells_detected"] = len({bc for barcodes in cell_barcodes.values() for bc in barcodes})

    return summary_info


def main(args, outs):
    with open(args.peaks, 'r') as infile:
        peak_regions = get_target_regions(infile)

    barcode_counts = Counter()
    targeted_counts = Counter()
    fragment_depth = 0
    for contig, start, stop, barcode, dup_count in parsed_fragments_from_contig(args.contig, args.fragments,
                                                                                args.fragments_index):
        barcode_counts[barcode] += 1
        fragment_depth += 1
        if fragment_overlaps_target(contig, start, stop, peak_regions):
            targeted_counts[barcode] += 1

    with open(outs.barcode_counts, 'w') as outfile:
        pickle.dump(barcode_counts, outfile)
    with open(outs.targeted_counts, 'w') as outfile:
        pickle.dump(targeted_counts, outfile)
    outs.fragment_depth = fragment_depth
