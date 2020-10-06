from __future__ import division

import scipy
import utils
from tenkit.stats import robust_divide
import numpy as np
from constants import NO_BARCODE

def add_bulk_targeting_metrics(summary_info, singlecell_df, species_list):
    """Take singlecell targeting data and calculate bulk targeting metrics from them.
    """
    for species in species_list:
        species_cell_mask = singlecell_df["is_%s_cell_barcode" % species] == 1
        key_suffix = "" if len(species_list) == 1 else "_{}".format(species)

        total = singlecell_df[species_cell_mask]["passed_filters"].sum()
        tss = singlecell_df[species_cell_mask]["TSS_fragments"].sum()
        dnase = singlecell_df[species_cell_mask]["DNase_sensitive_region_fragments"].sum()
        enhancer = singlecell_df[species_cell_mask]["enhancer_region_fragments"].sum()
        promoter = singlecell_df[species_cell_mask]["promoter_region_fragments"].sum()
        ontarget = singlecell_df[species_cell_mask]["on_target_fragments"].sum()
        blacklist = singlecell_df[species_cell_mask]["blacklist_region_fragments"].sum()
        peaks = singlecell_df[species_cell_mask]["peak_region_fragments"].sum()
        summary_info['frac_fragments_overlapping_targets{}'.format(key_suffix)] = robust_divide(ontarget, total)
        summary_info['frac_fragments_overlapping_tss{}'.format(key_suffix)] = robust_divide(tss, total)
        summary_info['frac_fragments_overlapping_dnase{}'.format(key_suffix)] = robust_divide(dnase, total)
        summary_info['frac_fragments_overlapping_enhancer{}'.format(key_suffix)] = robust_divide(enhancer, total)
        summary_info['frac_fragments_overlapping_promoter{}'.format(key_suffix)] = robust_divide(promoter, total)
        summary_info['frac_fragments_overlapping_blacklist{}'.format(key_suffix)] = robust_divide(blacklist, total)
        summary_info['frac_fragments_overlapping_peaks{}'.format(key_suffix)] = robust_divide(peaks, total)
    cell_mask = singlecell_df['cell_id'] != 'None'
    cut_frags_in_peaks = singlecell_df[cell_mask]["peak_region_cutsites"].sum()
    total = singlecell_df[cell_mask]["passed_filters"].sum()
    summary_info['frac_cut_fragments_in_peaks'] = robust_divide(cut_frags_in_peaks, 2 * total)

    return summary_info

def add_singlecell_sensitivity_metrics(summary_info, singlecell_df, species_list):
    """Get sensitivity metrics per cell"""

    for species in species_list:
        key_suffix = "" if len(species_list) == 1 else "_{}".format(species)

        cell_mask = singlecell_df["is_%s_cell_barcode" % species] == 1
        cell_fragments = singlecell_df['passed_filters'][cell_mask].values
        summary_info['mean_fragments_per_cell{}'.format(key_suffix)] = np.mean(cell_fragments)
        summary_info['median_fragments_per_cell{}'.format(key_suffix)] = np.median(cell_fragments)
        summary_info['min_fragments_per_cell{}'.format(key_suffix)] = np.min(cell_fragments) if cell_fragments.size else None
        summary_info['stdev_fragments_per_cell{}'.format(key_suffix)] = np.std(cell_fragments)
        cell_frags_overlapping_peaks = singlecell_df['peak_region_fragments{}'.format(key_suffix)][cell_mask].values
        summary_info['median_frags_overlapping_peaks_per_cell{}'.format(key_suffix)] = np.median(cell_frags_overlapping_peaks)

    # Background noise metric calculations
    valid_noncell_mask = (singlecell_df['cell_id'] == 'None') & (singlecell_df['barcode'] != NO_BARCODE)
    noncell_fragments = singlecell_df['passed_filters'][valid_noncell_mask].values
    summary_info['mean_fragments_per_noncell'] = np.mean(noncell_fragments)
    summary_info['median_fragments_per_noncell'] = np.median(noncell_fragments)

    return summary_info

def add_doublet_rate_metrics(summary_info, singlecell_df, species_list):
    """Infer doublet rate from observed doublets"""

    def infer_multiplets_from_observed(n_obs_multiplets, n_cells0, n_cells1):
        """Estimates the number of real multiplets based on the number observed from a barnyard (mixed species) experiment"""
        if n_cells0 == 0 or n_cells1 == 0 or n_obs_multiplets == 0:
            return 0

        # Prior probability of a doublet given counts for each cell type (ignore N_cells > 2)
        p_obs_multiplet = (2 * (n_cells0 / (n_cells0 + n_cells1)) * (n_cells1 / (n_cells0 + n_cells1)))

        # Brute force MLE of binomial n
        likelihood = scipy.stats.binom.pmf(n_obs_multiplets, xrange(0, n_cells0 + n_cells1), p_obs_multiplet)
        return np.argmax(likelihood)

    has_species_info = (species_list != [""])
    if not has_species_info or len(species_list) < 2:
        return summary_info

    counts = []
    cell_barcodes_dict = {}
    for species in species_list:
        species_cell_mask = singlecell_df["is_%s_cell_barcode" % species] == 1
        print singlecell_df['barcode'][species_cell_mask].values.tolist()
        cell_barcodes_dict[species] = singlecell_df['barcode'][species_cell_mask].values.tolist()
        counts.append(len(cell_barcodes_dict[species]))
    total_unique_cell_barcodes = {bc for barcodes in cell_barcodes_dict.values() for bc in barcodes}
    total_cell_barcodes = sum(counts)
    summary_info['cells_detected'] = len(total_unique_cell_barcodes)
    if len(species_list) > 1:
        observed_doublets = total_cell_barcodes - len(total_unique_cell_barcodes)
        observed_doublet_rate = robust_divide(observed_doublets, total_cell_barcodes)
        inferred_doublets = infer_multiplets_from_observed(observed_doublets, counts[0], counts[1])
        inferred_doublet_rate = robust_divide(inferred_doublets, total_cell_barcodes)
        summary_info['observed_doublets'] = observed_doublets
        summary_info['observed_doublet_rate'] = observed_doublet_rate
        summary_info['inferred_doublets'] = inferred_doublets
        summary_info['inferred_doublet_rate'] = inferred_doublet_rate

    return summary_info

def add_purity_metrics(summary_info, singlecell_df, species_list):
    """Detect purity of each barcode"""

    if len(species_list) == 2:
        spec_1_purity, spec_2_purity = utils.get_purity_info(singlecell_df, species_list)

        summary_info['median_purity_{}'.format(species_list[0])] = np.median(spec_1_purity)
        summary_info['median_purity_{}'.format(species_list[1])] = np.median(spec_2_purity)

    return summary_info
