"""
Tools for generating ATAC-specific plots from terminal pipeline outputs.
"""

from __future__ import division, print_function

import pandas as pd
from itertools import izip
from plots import downsample_scatterplot_by_density
from utils import get_purity_info
import numpy as np
from scipy import stats
from barcodes import load_barcode_whitelist
from tools.io import load_insert_sizes
from constants import NO_BARCODE, MINIMUM_COUNT, WHITELIST_CONTAM_RATE
from metrics import calculate_tss_score_and_profile, calculate_ctcf_score_and_profile

def generate_ctcf_enrichment_plot(ctcf_relpos, ctcf_file_path, debug):
    _, yvals, xvals = calculate_ctcf_score_and_profile(ctcf_relpos)

    data = {
        "layout": {
            "xaxis": {
                "type": "linear",
                "title": "Relative Position (bp from CTCF motif)",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "linear",
                "title": "Relative Enrichment",
                "showline": True,
                "zeroline": False,
            },
            "title": "Enrichment around CTCF motifs (normalized)",
        },
        "data": [
            {
                "name": "Cells",
                "x": list(xvals),
                "y": list(yvals),
                "type": "scatter",
                "mode": "lines",
            },
        ],
    }

    if debug:
        ctcf_ref = pd.read_csv(ctcf_file_path)
        xvals_ref = ctcf_ref['position']
        yvals_ref = ctcf_ref['cutsites']
        # => was already normalized.
        ref_dict = {"name": "Reference",
                    "x": list(xvals_ref),
                    "y": list(yvals_ref),
                    "type": "scatter",
                    "mode": "lines",
                    "marker": {"opacity": 0.7,
                               "color": "grey",
                               "size": 4}}
        data["data"].append(ref_dict)
    return data


def generate_tss_enrichment_plot(tss_relpos, tss_file_path, debug):
    _, yvals, xvals = calculate_tss_score_and_profile(tss_relpos)

    data = {
        "layout": {
            "xaxis": {
                "type": "linear",
                "title": "Relative Position (bp from TSS)",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "linear",
                "title": "Relative Enrichment",
                "showline": True,
                "zeroline": False,
            },
            "title": "Enrichment around TSS (normalized)",
        },
        "data": [
            {
                "name": "Cells",
                "x": list(xvals),
                "y": list(yvals),
                "type": "scatter",
                "mode": "lines",
            },
        ],
    }
    if debug:
        tss_ref = pd.read_csv(tss_file_path)
        xvals_ref = tss_ref['position']
        yvals_ref = tss_ref['cutsites']
        ref_dict = {"name": "Reference",
                    "x": list(xvals_ref),
                    "y": list(yvals_ref),
                    "type": "scatter",
                    "mode": "lines",
                    "marker": {"opacity": 0.7,
                               "color": "grey",
                               "size": 4}}
        data["data"].append(ref_dict)
    return data


def generate_targeting_plot(singlecell_df, species_list, target="TSS", excluded_barcodes=None):
    is_not_cell = (singlecell_df['cell_id'] == 'None') & (singlecell_df['barcode'] != NO_BARCODE)
    if len(species_list) == 1:
        labels = ["Non-cells", "Cells"]
        are_cells = [singlecell_df['is__cell_barcode'] == 1]
    else:
        labels = ["Non-cells"] + species_list
        are_cells = [singlecell_df['is_{}_cell_barcode'.format(species)] == 1
                     for species in species_list]
    masks = [is_not_cell] + are_cells
    total = singlecell_df["passed_filters"]
    target_map = {
        "TSS": "TSS_fragments",
        "Peaks": "peak_region_fragments",
        "DNAse": "DNase_sensitive_region_fragments",
        "Enhancer": "enhancer_region_fragments",
        "Promoter": "promoter_region_fragments",
        "On Target": "on_target_fragments",
    }

    # Specially show excluded barcode categories
    if excluded_barcodes is not None:
        exclusions_by_type = {}
        for species, barcode_data in excluded_barcodes.iteritems():
            for barcode in barcode_data:
                reason = barcode_data[barcode][0]
                if reason not in exclusions_by_type:
                    exclusions_by_type[reason] = set()
                exclusions_by_type[reason].add(barcode)
        for reason, barcode_set in exclusions_by_type.iteritems():
            labels.append(reason)
            exclude_mask = np.array([bc in barcode_set for bc in singlecell_df['barcode']])
            masks.append(exclude_mask)
            # Remove these from the "Non-cells" category as well
            masks[0][exclude_mask] = False

    no_dups_dfs = []
    for label, mask in zip(labels, masks):
        dups_df = pd.DataFrame({'total': total[mask], 'subtype': singlecell_df[target_map[target]][mask]})
        dups_df.drop_duplicates(inplace=True)
        no_dups_dfs.append(dups_df)
    if sum(map(len, no_dups_dfs)) > 2000 * len(species_list):
        downsampled_dfs = []
        for dataframe in no_dups_dfs:
            downsampled_dfs.append(downsample_scatterplot_by_density(dataframe, 2000, "total", "subtype"))
        no_dups_dfs = downsampled_dfs

    data = []
    for label, points_df, color in izip(labels, no_dups_dfs, ["red", "blue", "cyan", "magenta", "green", "orange"]):
        xvals = points_df['total']
        yvals = points_df['subtype'] / points_df['total']
        data.append({"name": label, "x": list(xvals), "y": list(yvals), "type": "scatter",
                     "mode": "markers",
                     "marker": {"opacity": 0.7,
                                "color": color,
                                "size": 4}})
    return {
        "layout": {
            "xaxis": {
                "type": "log",
                "title": "Fragments per Barcode",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "linear",
                "title": "Fraction Fragments Overlapping {}".format(target),
                "showline": True,
                "zeroline": False,
            },
            "title": "Singlecell Targeting ({})".format(target),
        },
        "data": data,
    }


def generate_clustering_plot(singlecell_df, tsne_results, cluster_labels, barcodes, species_list,
                             library_list=None, method='cluster'):
    def get_data_by_cluster(index, labels, name=None):
        if name is None:
            name = "Cluster {}".format(index)
        mask = labels == index
        name = "{} ({})".format(name, sum(mask))
        return {"name": name,
                "x": list(tsne_results[mask, 0]),
                "y": list(tsne_results[mask, 1]),
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.9, "size": 4},
                "text": "{} cells, {:.1%}".format(mask.sum(), mask.sum() / len(mask))}

    def get_data_by_depth():
        full_barcodes = singlecell_df['barcode'].values.tolist()
        bc_indices = np.array([full_barcodes.index(bc) for bc in barcodes], dtype=int)
        depths = singlecell_df["passed_filters"].values[bc_indices]
        return {
            "name": "Cells",
            "x": list(tsne_results[:, 0]),
            "y": list(tsne_results[:, 1]),
            "type": "scatter",
            "mode": "markers",
            "marker": {"opacity": 0.6,
                       "size": 6,
                       "cmin": np.log10(min(depths)),
                       "cmax": np.log10(max(depths)),
                       "color": list(np.log10(depths)),
                       "colorscale": "Viridis",
                       "colorbar": {"title": "log10 Fragments"}
                       },
            "text": ['{} fragments'.format(depth) for depth in depths],
        }

    def get_data_by_species():
        full_barcodes = singlecell_df['barcode'].values.tolist()
        bc_indices = np.array([full_barcodes.index(bc) for bc in barcodes], dtype=int)
        custom_labels = - np.ones(len(barcodes), dtype=int)
        mask_sp1 = singlecell_df["is_{}_cell_barcode".format(species_list[0])].values[bc_indices] == 1
        mask_sp2 = singlecell_df["is_{}_cell_barcode".format(species_list[1])].values[bc_indices] == 1
        custom_labels[mask_sp1] = 0
        custom_labels[mask_sp2] = 1
        custom_labels[mask_sp1 & mask_sp2] = 2
        names = species_list + ['Doublet']
        return [get_data_by_cluster(i, custom_labels, name) for i, name in enumerate(names)]

    def get_subdata(index, labels, text, name=None):
        if name is None:
            name = "{} {}".format(text, index)
        mask = labels == index
        name = "{} ({})".format(name, sum(mask))
        return {"name": name,
                "x": list(tsne_results[mask, 0]),
                "y": list(tsne_results[mask, 1]),
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.9, "size": 4},
                "text": "{} cells, {:.1%}".format(mask.sum(), mask.sum() / len(mask))}

    def get_data_by_gemgroup():
        gem_group_labels = np.array([int(bc.split("-")[1]) for bc in barcodes if "-" in bc])
        return [get_subdata(i + 1, gem_group_labels, "Library", name="Library {}".format(library_list[i]))
                for i in range(max(gem_group_labels))]

    if method == 'cluster':
        data = [get_data_by_cluster(i + 1, cluster_labels) for i in range(max(cluster_labels))]
    elif method == 'depth':
        data = [get_data_by_depth()]
    elif method == 'species':
        data = get_data_by_species()
    elif method == 'library':
        data = get_data_by_gemgroup()
    else:
        raise ValueError('Unexpected method provided: {}'.format(method))

    return {
        "layout": {
            "xaxis": {
                "type": "linear",
                "title": "tSNE axis 1",
                "showline": False,
                "zeroline": True,
            },
            "yaxis": {
                "type": "linear",
                "title": "tSNE axis 2",
                "showline": False,
                "zeroline": True,
            },
            "title": "Cell Clustering (By {}{})".format(method[0].upper(), method[1:]),
        },
        "data": data,
    }


def generate_knee_plot(singlecell_df, excluded_barcodes=None, species=""):
    # N.B. We call cells based on peak-targeted fragments, so that's what we show in the knee plot
    if species:
        fragment_counts = singlecell_df['peak_region_fragments_{}'.format(species)].values
    else:
        fragment_counts = singlecell_df['peak_region_fragments'].values

    if excluded_barcodes is None:
        valid_barcode_mask = singlecell_df['barcode'] != NO_BARCODE
    else:
        barcodes = singlecell_df['barcode'].values
        valid_barcode_mask = np.array([(bc not in excluded_barcodes[species]) and (bc != NO_BARCODE)
                                        for bc in barcodes])

    cell_mask = (singlecell_df['is_{}_cell_barcode'.format(species)] == 1) & valid_barcode_mask
    noncell_mask = (singlecell_df['is_{}_cell_barcode'.format(species)] == 0) & valid_barcode_mask
    data_subplots = []
    yvals = np.unique(fragment_counts[noncell_mask])
    yvals.sort()
    xshift = sum(cell_mask)
    xvals = [(fragment_counts[noncell_mask] >= count).sum() + xshift for count in yvals]
    data_subplots.append({
        "name": "Non-cells",
        "x": list(xvals),
        "y": list(yvals),
        "type": "scatter",
        "mode": "lines",
        "line": {"width": 4, "color": "blue"},
    })

    yvals = np.unique(fragment_counts[cell_mask])
    yvals.sort()
    xvals = [(fragment_counts[cell_mask] >= count).sum() for count in yvals]
    data_subplots.append({
        "name": "{} Cells".format(species),
        "x": list(xvals),
        "y": list(yvals),
        "type": "scatter",
        "mode": "lines",
        "line": {"width": 4, "color": "orange"},
    })

    return {
        "layout": {
            "xaxis": {
                "type": "log",
                "title": "Barcodes",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "log",
                "title": "{} Fragments Overlapping Peaks".format(species),
                "showline": True,
                "zeroline": False,
            },
            "title": "{} Cells".format(species),
        },
        "data": data_subplots,
    }


def generate_cellcalling_fragment_counts_plot(
        singlecell_df, cell_parameters, barcode_whitelist, excluded_barcodes=None, species=""):
    if species:
        fragment_counts = np.array(singlecell_df['passed_filters_{}'.format(species)].values)
    else:
        fragment_counts = np.array(singlecell_df['peak_region_fragments'].values)

    if excluded_barcodes is None:
        valid_barcode_mask = singlecell_df['barcode'] != NO_BARCODE
    else:
        barcodes = singlecell_df['barcode'].values
        valid_barcode_mask = np.array([(bc not in excluded_barcodes[species]) and (bc != NO_BARCODE) 
                                       for bc in barcodes], dtype=bool)

    threshold = cell_parameters[species]['cell_threshold']

    cell_mask = (fragment_counts >= threshold) & valid_barcode_mask
    noncell_mask = (fragment_counts < threshold) & valid_barcode_mask

    logbinmax = int(np.ceil(np.log10(fragment_counts.max())))
    xbins = list(np.hstack([np.arange(100), np.logspace(np.log10(100), logbinmax, 350)]))

    data_subplots = []
    for name, mask in zip(["Non-cells", "{} Cells".format(species)], [noncell_mask, cell_mask]):
        if mask.sum() > 0:
            counts, _ = np.histogram(fragment_counts[mask], xbins)
            data_subplots.append({
                "name": name,
                "x": xbins,
                "y": list(counts),
                "type": "scatter",
                "connectgaps": True,
                "fill": "tozeroy",
            })

    whitelist_length = len(load_barcode_whitelist(barcode_whitelist))
    fragment_depth = sum(singlecell_df['passed_filters'].values)
    count_shift = max(MINIMUM_COUNT, int(fragment_depth * WHITELIST_CONTAM_RATE / whitelist_length))

    def get_fitted_counts(barcode_total, bins, species, parameters):
        max_count = max(bins)
        count_values = np.arange(max_count + 1)

        frac_noise = parameters[species]['fraction_noise']
        mean_noise = parameters[species]['noise_mean']
        mean_signal = parameters[species]['signal_mean']
        dispersion_noise = parameters[species]['noise_dispersion']
        dispersion_signal = parameters[species]['signal_dispersion']

        estimated_noise_counts = stats.nbinom.pmf(count_values,
                                                  1 / dispersion_noise,
                                                  1 / (1 + dispersion_noise * mean_noise))
        estimated_signal_counts = stats.nbinom.pmf(count_values,
                                                   1 / dispersion_signal,
                                                   1 / (1 + dispersion_signal * mean_signal))
        estimated_noise_counts *= frac_noise * barcode_total
        estimated_signal_counts *= (1 - frac_noise) * barcode_total
        noise_bin_counts = np.array([estimated_noise_counts[(count_values >= lower) & (count_values < upper)].sum()
                                     for lower, upper in zip(bins[:-1], bins[1:])])
        signal_bin_counts = np.array([estimated_signal_counts[(count_values >= lower) & (count_values < upper)].sum()
                                     for lower, upper in zip(bins[:-1], bins[1:])])
        noise_bin_counts[noise_bin_counts < 1.0] = 0.0
        signal_bin_counts[signal_bin_counts < 1.0] = 0.0
        return bins[:-1], noise_bin_counts, signal_bin_counts

    xvals, noise, signal = get_fitted_counts((fragment_counts >= count_shift).sum(), xbins, species, cell_parameters)
    data_subplots.append({
        "name": "Noise fit",
        "x": list(xvals),
        "y": list(noise),
        "type": "scatter",
        "mode": "lines",
        "line": {"color": "grey",
                 "width": 1},
    })
    data_subplots.append({
        "name": "Signal fit",
        "x": list(xvals),
        "y": list(signal),
        "type": "scatter",
        "mode": "lines",
        "line": {"color": "black",
                 "width": 1},
    })
    data_subplots.append({
        "name": "Joint fit",
        "x": list(xvals),
        "y": list(signal + noise),
        "type": "scatter",
        "mode": "lines",
        "line": {"color": "red",
                 "width": 1},
    })

    return {
        "layout": {
            "xaxis": {
                "type": "log",
                "title": "{} Fragments Per Barcode".format(species),
            },
            "yaxis": {
                "type": "log",
                "title": "Barcodes",
            },
            "title": "{} Fragment Distribution".format(species),
        },
        "data": data_subplots,
    }


def generate_fragment_counts_plot(singlecell_df, species=""):
    if species:
        fragment_counts = singlecell_df['passed_filters_{}'.format(species)].values
    else:
        fragment_counts = singlecell_df['passed_filters'].values
    cell_mask = (singlecell_df['is_{}_cell_barcode'.format(species)] == 1) & (singlecell_df['barcode'] != NO_BARCODE)
    noncell_mask = (singlecell_df['is_{}_cell_barcode'.format(species)] == 0) & (singlecell_df['barcode'] != NO_BARCODE)

    logbinmax = np.ceil(np.log10(fragment_counts.max()))
    xbins = list(np.hstack([np.arange(100), np.logspace(np.log10(100), logbinmax, 350)]))

    data_subplots = []
    for name, mask, color in zip(["Non-cells", "{} Cells".format(species)],
                                 [noncell_mask, cell_mask],
                                 ["blue", "orange"]):
        # protect against really low depth samples, mixed up ref, severe under cell calling
        if mask.sum() > 0 and logbinmax > np.log10(100):
            counts, _ = np.histogram(fragment_counts[mask], xbins)
            data_subplots.append({
                "name": name,
                "x": xbins,
                "y": list(counts),
                "type": "scatter",
                "connectgaps": True,
                "fill": "tozeroy",
                "line": {"color": color},
            })
    return {
        "layout": {
            "xaxis": {
                "type": "log",
                "title": "{} Fragments Per Barcode".format(species),
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "log",
                "title": "Barcodes",
                "showline": True,
                "zeroline": False,
            },
            "title": "{} Fragment Distribution".format(species),
        },
        "data": data_subplots,
    }


def generate_wasted_data_plot(singlecell_df, title="Wasted Fragments Flow Diagram"):
    total_fragments = singlecell_df['total'].sum()
    is_cell_mask = singlecell_df['cell_id'] != 'None'
    valid_noncell_mask = (singlecell_df['cell_id'] == 'None') & (singlecell_df['barcode'] != NO_BARCODE)
    invalid_bc_mask = singlecell_df['barcode'] == NO_BARCODE
    invalid_barcodes = singlecell_df['total'][invalid_bc_mask].sum()
    noncell_barcodes = singlecell_df['total'][valid_noncell_mask].sum()
    cell_barcodes = singlecell_df['total'][is_cell_mask].sum()

    # Split cell barcode signal up into various wasted types
    cell_mito = singlecell_df[is_cell_mask]['mitochondrial'].sum()
    cell_unmapped = singlecell_df[is_cell_mask]['unmapped'].sum()
    cell_lowmapq = singlecell_df[is_cell_mask]['lowmapq'].sum()
    cell_dups = singlecell_df[is_cell_mask]['duplicate'].sum()
    cell_chimeric = singlecell_df[is_cell_mask]['chimeric'].sum()
    cell_passed = singlecell_df[is_cell_mask]['passed_filters'].sum()

    # Split up non-cell barcode signal as well
    noncell_mito = singlecell_df[valid_noncell_mask]['mitochondrial'].sum()
    noncell_unmapped = singlecell_df[valid_noncell_mask]['unmapped'].sum()
    noncell_lowmapq = singlecell_df[valid_noncell_mask]['lowmapq'].sum()
    noncell_dups = singlecell_df[valid_noncell_mask]['duplicate'].sum()
    noncell_chimeric = singlecell_df[valid_noncell_mask]['chimeric'].sum()
    noncell_passed = singlecell_df[valid_noncell_mask]['passed_filters'].sum()

    return {
        "layout": {
            "title": title,
            "height": 800,
        },
        "data": [{
            "type": "sankey",
            "orientation": "h",
            "arrangement": "snap",
            "node": {
                "color": [
                    "blue",
                    "orange",
                    "red",
                    "green",
                    "pink",
                    "pink",
                    "pink",
                    "pink",
                    "pink",
                    "magenta",
                    "cyan",
                    "cyan",
                    "cyan",
                    "cyan",
                    "cyan",
                    "green",
                ],
                "pad": 30,
                "thickness": 80,
                "line": {"width": 0},
                "label": [
                    "Total Fragments",
                    "Invalid Barcodes ({:.1%})".format(invalid_barcodes / total_fragments),
                    "Non-cell Barcodes ({:.1%})".format(noncell_barcodes / total_fragments),
                    "Cell Barcodes ({:.1%})".format(cell_barcodes / total_fragments),
                    "Mitochondrial Fragments ({:.1%})".format(noncell_mito / total_fragments),
                    "Unmapped Fragments ({:.1%})".format(noncell_unmapped / total_fragments),
                    "Low MapQ Fragments ({:.1%})".format(noncell_lowmapq / total_fragments),
                    "Duplicate Fragments ({:.1%})".format(noncell_dups / total_fragments),
                    "Chimeric Fragments ({:.1%})".format(noncell_chimeric / total_fragments),
                    "High Quality Fragments ({:.1%})".format(noncell_passed / total_fragments),
                    "Mitochondrial Fragments ({:.1%})".format(cell_mito / total_fragments),
                    "Unmapped Fragments ({:.1%})".format(cell_unmapped / total_fragments),
                    "Low MapQ Fragments ({:.1%})".format(cell_lowmapq / total_fragments),
                    "Duplicate Fragments ({:.1%})".format(cell_dups / total_fragments),
                    "Chimeric Fragments ({:.1%})".format(cell_chimeric / total_fragments),
                    "High Quality Fragments ({:.1%})".format(cell_passed / total_fragments),
                ],
            },
            "link": {
                "source": [0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3],
                "target": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
                "value": [invalid_barcodes, noncell_barcodes, cell_barcodes,
                          noncell_mito, noncell_unmapped, noncell_lowmapq,
                          noncell_dups, noncell_chimeric, noncell_passed,
                          cell_mito, cell_unmapped, cell_lowmapq, cell_dups, cell_chimeric, cell_passed]
            }
        }]
    }


def generate_insertsize_plot(singlecell_df, insertsize_fn, scale='linear'):
    xvals = np.arange(25, 800)
    cell_barcodes = set(singlecell_df['barcode'][singlecell_df['cell_id'] != 'None'].values)
    yvals = load_insert_sizes(insertsize_fn, colrange=xvals, barcode_selection=cell_barcodes, report_sum=True)
    return {
        "layout": {
            "xaxis": {
                "type": "linear",
                "title": "Insert Size",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": scale,
                "title": "Fragment Count ({} scale)".format(scale),
                "showline": True,
                "zeroline": False,
            },
            "title": "Insert Size Distribution",
        },
        "data": [
            {
                "name": "Cells",
                "x": list(xvals),
                "y": list(yvals),
                "type": "scatter",
                "mode": "lines",
                "line": {"color": "blue"},
            },
        ],
    }


def generate_purity_plot(singlecell_df, species_list):
    """Takes a dataframe of the singlecell.csv output of the ATAC pipeline and a list of
    species and produces an output dictionary of data for a purity plot."""
    if len(species_list) < 2:
        raise ValueError('Cannot generate a purity plot on single-species data')

    spec_1_purity, spec_2_purity = get_purity_info(singlecell_df, species_list)
    purity_data = []
    for species, purity in zip(species_list, [spec_1_purity, spec_2_purity]):
        yvals, xvals = np.histogram(purity, bins=25, range=(0.5, 1))
        xvals = (xvals[:-1] + xvals[1:]) / 2
        purity_data.append({
            "name": "{} Cells".format(species),
            "x": list(xvals),
            "y": list(yvals),
            "type": "scatter",
            "fill": "tozeroy",
        })

    return {
        "layout": {
            "xaxis": {
                "type": "linear",
                "title": "Fraction of Fragments in Primary Species",
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "linear",
                "title": "# of Cell Barcodes",
                "showline": True,
                "zeroline": False,
            },
            "title": "Purity (Cell Barcodes)",
        },
        "data": purity_data
    }


def generate_barnyard_plot(singlecell_df, species_list, downsample=True):
    """Takes a dataframe of the singlecell.csv output of the ATAC pipeline and a list of
    species and produces an output dictionary of barnyard plot data.
    """
    if len(species_list) < 2:
        raise ValueError('Cannot generate a barnyard plot on single-species data')

    xdata = singlecell_df["passed_filters_{}".format(species_list[0])].values
    ydata = singlecell_df["passed_filters_{}".format(species_list[1])].values
    is_not_cell = singlecell_df['cell_id'] == 'None'
    is_spec_1 = singlecell_df['is_{}_cell_barcode'.format(species_list[0])] == 1
    is_spec_2 = singlecell_df['is_{}_cell_barcode'.format(species_list[1])] == 1
    multiplets = is_spec_1 & is_spec_2
    is_spec_1[multiplets] = False
    is_spec_2[multiplets] = False

    masks = [is_not_cell, multiplets, is_spec_1, is_spec_2]
    labels = ['noncells', 'multiplets', species_list[0], species_list[1]]

    # remove duplicate points.
    no_dups_dfs = []
    for mask in masks:
        dups_df = pd.DataFrame({'x': xdata[mask], 'y': ydata[mask]})
        dups_df.drop_duplicates(inplace=True)
        no_dups_dfs.append(dups_df)

    # if more than 2k*(numb of groups) points in all groups, and downsample enabled,
    # downsample to 2k*(numb of groups) points by sampling 2000 points in each group.
    # if downsample not enabled, only sample 2000 points in the noncells group.
    if sum(map(len, no_dups_dfs)) > 2000 * (1 + len(species_list)):
        for label, df in izip(labels, no_dups_dfs):
            df["label"] = label
        ds_noncells = downsample_scatterplot_by_density(no_dups_dfs[0], 2000, 'x', 'y')

        # multiplets will not be downsampled.
        multiplets_df = no_dups_dfs[1]

        if downsample:
            full_cells_df = pd.concat(no_dups_dfs[2:])
            full_cells_df.reset_index(drop=True, inplace=True)
            full_cells_df_ds = downsample_scatterplot_by_density(full_cells_df, 2000 * len(species_list), 'x', 'y')
            ds_cells = [full_cells_df_ds.loc[full_cells_df_ds['label'] == label, ['x', 'y']]
                        for label in labels[2:]]
        else:
            ds_cells = no_dups_dfs[2:]

        no_dups_dfs = [ds_noncells[['x', 'y']]] + [multiplets_df[['x', 'y']]] + ds_cells

    xvals = {'noncells': list(no_dups_dfs[0]["x"]),
             'multiplets': list(no_dups_dfs[1]["x"]),
             species_list[0]: list(no_dups_dfs[2]["x"]),
             species_list[1]: list(no_dups_dfs[3]["x"]),
             }
    yvals = {'noncells': list(no_dups_dfs[0]["y"]),
             'multiplets': list(no_dups_dfs[1]["y"]),
             species_list[0]: list(no_dups_dfs[2]["y"]),
             species_list[1]: list(no_dups_dfs[3]["y"]),
             }

    return {
        "layout": {
            "xaxis": {
                "type": "log",
                "title": "{} Fragments Per Cell".format(species_list[0]),
                "showline": True,
                "zeroline": False,
            },
            "yaxis": {
                "type": "log",
                "title": "{} Fragments Per Cell".format(species_list[1]),
                "showline": True,
                "zeroline": False,
            },
            "title": "Barnyard",
        },
        "data": [
            {
                "name": "Non-cells",
                "x": xvals['noncells'],
                "y": yvals['noncells'],
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.4,
                           "color": "red",
                           "size": 3},
            },
            {
                "name": "Multiplets",
                "x": xvals['multiplets'],
                "y": yvals['multiplets'],
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.7,
                           "color": "green",
                           "size": 6},
            },
            {
                "name": "{} Cells".format(species_list[0]),
                "x": xvals[species_list[0]],
                "y": yvals[species_list[0]],
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.5,
                           "color": "blue",
                           "size": 4},
            },
            {
                "name": "{} Cells".format(species_list[1]),
                "x": xvals[species_list[1]],
                "y": yvals[species_list[1]],
                "type": "scatter",
                "mode": "markers",
                "marker": {"opacity": 0.5,
                           "color": "cyan",
                           "size": 4},
            },
        ]

    }
