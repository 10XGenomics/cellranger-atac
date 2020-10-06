"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Produces the final customer-facing web_summary.
"""

from __future__ import division

from websummary import summarize

import json
import pandas as pd
import martian
import numpy as np
import os
import h5py

from tools import ReferenceManager
from metrics import MetricAnnotations

import cellranger.matrix as cr_matrix
import cellranger.analysis.constants as cr_analysis_constants
import cellranger.analysis.graphclust as cr_graphclust

from constants import RPC_30K, RPC_50K, RPC_10K

from plots.plotly_tools import PLOT_CONFIG_KWARGS
from plots.atac import (generate_barnyard_plot, generate_purity_plot, generate_insertsize_plot,
                        generate_wasted_data_plot, generate_fragment_counts_plot, generate_knee_plot,
                        generate_clustering_plot, generate_targeting_plot, generate_tss_enrichment_plot,
                        generate_ctcf_enrichment_plot)

__MRO__ = """
stage CREATE_WEBSUMMARY(
    in  string reference_path,
    in  string barcode_whitelist,
    in  json   summary_results,
    in  json   bulk_complexity,
    in  json   singlecell_complexity,
    in  string sample_id,
    in  string sample_desc,
    in  map[]  sample_def,
    in  bool   debug,
    in  csv    singlecell,
    in  csv    insert_sizes,
    in  csv    tss_relpos,
    in  csv    ctcf_relpos,
    in  h5     filtered_peak_bc_matrix,
    in  h5     analysis,
    in  json   excluded_barcodes,
    out html   web_summary,
    src py     "stages/reporter/create_websummary",
) using (
    mem_gb = 16,
)
"""

DOWNSAMPLE_BARNYARD = True
DOWNSAMPLE_TARGETING = True


def get_hero_metric_data(metadata, summary_data, species_list, debug):
    """Adds a few metrics of primary importance into the header of the web summary:
    - Annotated cell counts (split by species if barnyard)
    - Median fragments per cell (split)
    - Non-dup wasted data
    Alerts raised depend on debug status.
    """
    data = {}
    if summary_data is None:
        return data

    for key in ['annotated_cells',
                'median_fragments_per_cell',
                'frac_fragments_overlapping_targets',
                'frac_cut_fragments_in_peaks']:
        metrics = metadata.gen_metric_list(summary_data, [key], species_list=species_list, debug=debug)
        for metric in metrics:
            key = metric.key
            if len(species_list) > 1:
                for i, species in enumerate(species_list):
                    key = key.replace(species, 'spec{}'.format(i + 1))
            data[key] = metric.gen_metric_dict()

            # low number of annotated_cells is a special case for metric coloring.
            # greater than case is handled automatically with the metrics.csv ranges;
            # less than case is handled here.
            # the 'translate' call removes possible thousand comma separators.
            if 'annotated_cells' in key and int(data[key]['metric'].translate(None, ',')) < 500:
                data[key]['threshold'] = "warn"

    # Alerts.
    alarm_keys = [
        'annotated_cells',
        'median_fragments_per_cell',
        'frac_fragments_overlapping_targets',
        'frac_cut_fragments_in_peaks'
    ]
    alarms = metadata.gen_metric_list(summary_data, alarm_keys, species_list=species_list, debug=debug)
    new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    if new_alarms:
        data['alarms'] = new_alarms

    # low number of annotated_cells is a special alert case (see comment above).

    for species in species_list:
        suffix = "_{}".format(species) if len(species_list) > 1 else ""
        suffix_withpar = " ({})".format(species) if len(species_list) > 1 else ""
        sufkey = "annotated_cells{}".format(suffix)
        if sufkey in summary_data:
            if summary_data[sufkey] < 500:
                # create alarm
                alarm_dict = {
                    "raw_value": summary_data["annotated_cells{}".format(suffix)],
                    "formatted_value": '{:,.0f}'.format(summary_data["annotated_cells{}".format(suffix)]),
                    "raised": True,
                    "parent": "annotated_cells{}".format(suffix),
                    "title": "Estimated number of cells is low{}".format(suffix_withpar),
                    "message": "Number of cells detected is expected to be higher than 500.  This usually \
                                    indicates poor cell, library, or sequencing quality.",
                    "level": "WARN",
                    "test": "",
                    "id": "annotated_cells{}".format(suffix)
                }
                if 'alarms' in data:
                    data['alarms'].append(alarm_dict)
                else:
                    data['alarms'] = [alarm_dict]

    return data


def get_barnyard_data(metadata, summary_data, species_list, singlecell_df, debug, downsample):
    """Adds a barnyard plot for multi-species samples.  Alerts raised depend on debug status."""
    data = {}
    if summary_data is None:
        return data

    metric_keys = ['observed_doublet_rate',
                   'inferred_doublet_rate',
                   'median_purity']
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    rows = [[metric.name, metric.value_string] for metric in metrics]
    data['barnyard_metrics_table'] = {
        "rows": rows
    }

    helptext_metrics = metadata.gen_metric_helptext(metric_keys)
    helptext_plots = [
        ["Plots", ["(left) Barnyard scatter plot, where each dot represents a barcode and its coordinates indicate "
                   "number of fragments, from each species, assigned to the barcode. Groups are estimated "
                   "computationally.",
                   "(right) Histograms of per barcode purities estimated for barcodes in each species."]]
    ]
    helptext = helptext_metrics + helptext_plots
    data['barnyard_helptext'] = {'title': 'Barnyard', 'data': helptext}

    if singlecell_df is None:
        return data

    data['plot_barnyard'] = generate_barnyard_plot(singlecell_df, species_list, downsample)
    data['plot_purity'] = generate_purity_plot(singlecell_df, species_list)

    return data


def get_cell_metrics_data(metadata, summary_data, species_list, singlecell_df, excluded_barcode_fn, debug):
    """Adds a section with per-cell sensitivity metrics, and fragments per cell histograms.
    """
    data = {}
    if summary_data is None:
        return data

    metric_keys = ['annotated_cells',
                   'cell_threshold',
                   'median_fragments_per_cell',
                   'median_fragments_per_noncell']
    if debug:
        metric_keys.append('estimated_gelbead_doublet_rate')
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    data_key = 'cell_metrics_table'
    data[data_key] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    helptext_metrics = metadata.gen_metric_helptext(metric_keys)
    helptext_plots = [
        ["Plots", ["(left) Knee plot of number of fragments overlapping peaks for all the barcodes in the library. "
                   "This number is used to call cells.",
                   "(right) Histograms of number of fragments per cell barcode for non-cells and cells."]
         ]
    ]
    helptext = helptext_metrics + helptext_plots
    data['cell_helptext'] = {'title': 'Cells', 'data': helptext}

    if singlecell_df is None:
        return data

    if excluded_barcode_fn is not None:
        with open(excluded_barcode_fn, 'r') as infile:
            excluded_barcodes = json.load(infile)
    else:
        excluded_barcodes = None

    for i, species in enumerate(species_list):
        key_suffix = '_spec{}'.format(i + 1) if species else ''
        data['fragments_per_cell_plot{}'.format(key_suffix)] = generate_fragment_counts_plot(singlecell_df, species)
        data['barcode_knee_plot{}'.format(key_suffix)] = generate_knee_plot(singlecell_df, excluded_barcodes, species)

    return data


def get_wasted_data(metadata, summary_data, singlecell_df, species_list, debug):
    """Add metrics on wasted data fractions and ratios along with a
    Sankey flow diagram of fragments.  Alerts raised depend on debug status.
    """
    data = {}
    if summary_data is None:
        return data

    metric_keys = ['num_fragments',
                   'total_usable_fragments',
                   'frac_valid_barcode',
                   'frac_valid_noncell',
                   'frac_waste_overall_nondup',
                   'waste_ratio_mito_cells']
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    data['wasted_data_table'] = {
        "rows": [[metric.name, metric.value_string] for metric in metrics]
    }

    helptext_metrics = metadata.gen_metric_helptext(metric_keys)
    helptext_plots = [
        ["Plots",
         ["Sankey flow diagram showing the final label assigned to each of the read pairs sequenced in the library."]],
    ]
    helptext = helptext_metrics + helptext_plots
    data['wasted_helptext'] = {'title': 'Wasted Data', 'data': helptext}
    if singlecell_df is None:
        return data

    data["sankey_fragments_plot"] = generate_wasted_data_plot(singlecell_df)

    # Alerts.
    alarm_keys = ['frac_valid_barcode']
    alarms = metadata.gen_metric_list(summary_data, alarm_keys, species_list, debug=debug)
    new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    if new_alarms:
        data['alarms'] = new_alarms

    return data


def get_pipeline_info(args, reference, debug):
    """Generates a table of general pipeline information.
    """
    metadata = reference.metadata

    def get_fastq_paths(sample_def):
        if sample_def is None:
            return ""
        else:
            paths = [x["read_path"] for x in sample_def]
            return "\n".join(paths)

    rows = [
        ['Sample ID', args.sample_id],
        ['Sample description', args.sample_desc],
        ['FASTQ path', get_fastq_paths(args.sample_def)],
        ['Pipeline version', martian.get_pipelines_version()],
        ['Reference path', args.reference_path],
    ]

    if metadata:
        rows.extend([
            ['Organism', metadata.get('organism')],
            ['Assembly', metadata.get('assembly')],
            ['Annotation', metadata.get('annotation')],
        ])

    if debug:
        rows.append(['Barcode Whitelist', args.barcode_whitelist])

    data = {'pipeline_info_table': {'rows': rows}}
    data['pipeline_helptext'] = {'title': 'Sample', 'data': []}
    return data


def get_sequencing_info(metadata, summary_data, species_list, debug):
    """ Alerts raised depend on debug status.
    """
    data = {}
    if summary_data is None:
        return data

    metric_keys = [
        'num_fragments',
        'frac_valid_barcode',
        'r1_q30_bases_fract',
        'r2_q30_bases_fract',
        'bc_q30_bases_fract',
        'si_q30_bases_fract',
    ]
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    data['sequencing_info_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    # Alerts.
    # 'frac_valid_barcode' was already taken care of in get_wasted_data, don't need to redo it here.
    alarm_keys = [
        'r1_q30_bases_fract',
        'r2_q30_bases_fract',
        'bc_q30_bases_fract'
    ]
    alarms = metadata.gen_metric_list(summary_data, alarm_keys, species_list, debug=debug)
    new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    if new_alarms:
        data['alarms'] = new_alarms

    helptext = metadata.gen_metric_helptext(metric_keys)
    data['sequencing_helptext'] = {'title': 'Sequencing', 'data': helptext}

    return data


def get_insertsize_data(metadata, summary_data, singlecell_df, insertsize_fn, species_list, debug):
    data = {}
    if summary_data is None:
        return data

    if debug:
        metric_keys = [
            'insert_twist_period',
            'insert_nucleosome_period',
            'frac_fragments_nfr',
            'frac_fragments_nuc',
        ]
    else:
        metric_keys = [
            'frac_fragments_nfr',
            'frac_fragments_nuc',
        ]

    metrics = metadata.gen_metric_list(summary_data, metric_keys)
    data['insert_size_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    helptext_metrics = metadata.gen_metric_helptext(metric_keys)
    if debug:
        helptext_plots = [["Plots",
                           ["(left) Insert size distribution in log scale.",
                            "(right) Insert size distribution in linear scale."]]]
    else:
        helptext_plots = [["Plots",
                           ["Insert size distribution in linear scale."]]]

    helptext = helptext_metrics + helptext_plots
    data['insertsize_helptext'] = {'title': 'Insert Sizes', 'data': helptext}

    if singlecell_df is None or insertsize_fn is None:
        return data

    if debug:
        data['insert_size_plot'] = generate_insertsize_plot(singlecell_df, insertsize_fn, 'log')
    data['non_log_insert_size_plot'] = generate_insertsize_plot(singlecell_df, insertsize_fn, 'linear')

    return data


def get_targeting_data(metadata, summary_data, species_list, singlecell_df, tss_relpos, ctcf_relpos, debug, downsample):
    data = {}
    if summary_data is None:
        return data

    metric_keys = [
        'tss_enrichment_score',
        'frac_fragments_overlapping_tss',
        'frac_fragments_overlapping_peaks',
        'frac_cut_fragments_in_peaks',
        'frac_fragments_overlapping_targets',
    ]
    if debug:
        metric_keys += [
            'frac_fragments_overlapping_dnase',
            'frac_fragments_overlapping_enhancer',
            'frac_fragments_overlapping_promoter',
            'frac_fragments_overlapping_blacklist',
            'ctcf_enrichment_score'
        ]
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list=species_list, debug=debug)
    data['targeting_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}
    helptext_targeting = metadata.gen_metric_helptext(metric_keys)

    metric_keys = [
        'frac_mapped_confidently',
        'frac_waste_unmapped',
        'frac_waste_mitochondrial'
    ]
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    data['wasted_targeting_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}
    helptext_wasted = metadata.gen_metric_helptext(metric_keys)

    helptext_plots = [
        ["Plots",
         ["(left) TSS profile, as described above.",
          "(right) Targeting scatter plot. Each dot represents a barcode. Horizontal axis is the barcode's number of fragments, vertical axis is the percentage of those fragments that overlap peaks. Non-cell and cell groups are represented with different colors."]],
    ] if not debug else [
        ["Plots",
         ["(1st row, left) TSS profile, as described above.",
          "(1st row, right) CTCF profile, as described above.",
          "(all other rows) Targeting scatter plots. Each dot represents a barcode. Horizontal axis is the barcode's number of fragments, vertical axis is the percentage of those fragments that overlap different sets of functional genomic regions. Each panel corresponds to a different kind of region. Non-cell and cell groups are represented with different colors."]],
    ]

    helptext = helptext_targeting + helptext_wasted + helptext_plots

    data['targeting_helptext'] = {'title': 'Targeting', 'data': helptext}
    if singlecell_df is None:
        return data

    if tss_relpos is not None:
        tss_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tss_ref.csv')
        data['tss_enrichment_plot'] = generate_tss_enrichment_plot(tss_relpos, tss_file_path, debug)

    if ctcf_relpos is not None and debug:
        ctcf_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ctcf_ref.csv')
        data['ctcf_enrichment_plot'] = generate_ctcf_enrichment_plot(ctcf_relpos, ctcf_file_path, debug)

    data['peak_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "Peaks")
    if debug:
        data['singlecell_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "On Target")
        data['tss_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "TSS")
        data['dnase_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "DNAse")
        data['enhancer_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "Enhancer")
        data['promoter_targeting_plot'] = generate_targeting_plot(singlecell_df, species_list, "Promoter")

    return data


def get_peakcalling_data(metadata, summary_data, species_list, debug):
    data = {}
    if summary_data is None:
        return data

    metric_keys = {
        'total_peaks_detected',
        'frac_cut_fragments_in_peaks',
    }
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)

    data['peak_calling_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}
    helptext = metadata.gen_metric_helptext(metric_keys)
    data['peakcalling_helptext'] = {'title': 'Peak Calling', 'data': helptext}

    return data


def get_complexity_data(metadata, summary_data, bulk_complexity, singlecell_complexity, species_list, debug):
    data = {}
    if summary_data is None:
        return data

    metric_keys = {
        'bulk_total_library_complexity',
        'bulk_estimated_saturation',
        'frac_waste_duplicate',
    }
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    data['bulk_complexity_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}

    helptext_metrics = metadata.gen_metric_helptext(metric_keys)
    helptext_plots = [
        ["Plots",
         ["(top, left) Observed and fitted bulk library complexity as a function of downsampling rate.",
          "(top, right) Observed per cell complexity as a function of downsampling rate in mean reads per cell.",
          "(bottom) Observed per cell complexity as a function of downsampling rate in mean reads per cell, but only sampling high quality read pairs."]],
    ] if debug else [
        ["Plots",
         ["Observed per cell complexity as a function of downsampling rate in mean reads per cell."]],
    ]
    helptext = helptext_metrics + helptext_plots
    data['complexity_helptext'] = {'title': 'Library Complexity', 'data': helptext}

    # TODO: move complexity plots into atac library code
    if debug and bulk_complexity is not None:
        with open(bulk_complexity, 'r') as infile:
            complexity_data = json.load(infile)
        data['bulk_complexity_plot'] = {
            "layout": {
                "xaxis": {
                    "type": "linear",
                    "title": "Total Fragments",
                    "showline": True,
                    "zeroline": False,
                },
                "yaxis": {
                    "type": "linear",
                    "title": "Unique Fragments",
                    "showline": True,
                    "zeroline": False,
                },
                "title": "Library Complexity",
            },
            "data": [
                {
                    "name": "Observed",
                    "x": list(complexity_data['total']),
                    "y": list(complexity_data['unique']),
                    "type": "scatter",
                    "mode": "lines",
                },
            ],
        }
        # plot the estimate if it was successfully computed
        if 'bulk_total_library_complexity' in summary_data.keys():
            data['bulk_complexity_plot']['data'] += [
                {
                    "name": "Estimated Complexity",
                    "x": [0, max(complexity_data['plot_total'])],
                    "y": [summary_data['bulk_total_library_complexity'], summary_data['bulk_total_library_complexity']],
                    "type": "scatter",
                    "mode": "lines",
                },
                {
                    "name": "Fitted",
                    "x": list(complexity_data['plot_total']),
                    "y": list(complexity_data['plot_estimated_unique']),
                    "type": "scatter",
                    "mode": "lines",
                },
            ]

    if singlecell_complexity is not None:
        with open(singlecell_complexity, 'r') as infile:
            complexity_data = json.load(infile)

        if debug:
            data['singlecell_complexity_plot'] = {
                "layout": {
                    "xaxis": {
                        "type": "linear",
                        "title": "Mean High-Quality Reads Per Cell",
                        "showline": True,
                        "zeroline": False,
                    },
                    "yaxis": {
                        "type": "linear",
                        "title": "Median Unique Fragments Per Cell",
                        "showline": True,
                        "zeroline": False,
                    },
                    "title": "Per-Cell Library Complexity",
                },
                "data": [
                    {
                        "name": "Observed",
                        "x": [c * 2 for c in complexity_data['total']],
                        "y": list(complexity_data['unique']),
                        "type": "scatter",
                        "mode": "lines",
                    },
                ],
            }
            # plot the estimate if it was successfully computed
            if 'median_per_cell_total_library_complexity' in summary_data.keys():
                data['singlecell_complexity_plot']['data'] += [{
                    "name": "Estimated Complexity",
                    "x": [0, max(complexity_data['total']) * 2],
                    "y": [summary_data['median_per_cell_total_library_complexity'], summary_data['median_per_cell_total_library_complexity']],
                    "type": "scatter",
                    "mode": "lines",
                }]

        # plot frags/cell vs total library depth
        data['assay_efficiency_plot'] = {
            "layout": {
                "xaxis": {
                    "type": "linear",
                    "title": "Mean Reads Per Cell",
                    "showline": True,
                    "zeroline": False,
                },
                "yaxis": {
                    "type": "linear",
                    "title": "Median Unique Fragments Per Cell",
                    "showline": True,
                    "zeroline": False,
                },
                "title": "Per-Cell Library Complexity At Read Depth",
            },
            "data": [
                {
                    "name": "Observed",
                    "x": list(complexity_data['total_depth']),
                    "y": list(complexity_data['unique']),
                    "type": "scatter",
                    "mode": "lines",
                },
            ],
        }

    return data


def get_clustering_plots(analysis, filtered_peak_bc_matrix, species_list, singlecell_df, is_barnyard):
    data = {}
    if filtered_peak_bc_matrix is None or analysis is None:
        return data
    tsne_results = np.array(h5py.File(analysis, 'r')[cr_analysis_constants.ANALYSIS_H5_TSNE_GROUP]
                            ['_2']['transformed_tsne_matrix'])
    # N.B.: these cluster labels are 1-indexed
    cluster_labels = cr_graphclust.load_graphclust_from_h5(analysis).clusters
    barcodes = cr_matrix.CountMatrix.load_h5_file(filtered_peak_bc_matrix).bcs
    data['graph_cluster_plot'] = generate_clustering_plot(singlecell_df, tsne_results, cluster_labels, barcodes,
                                                          species_list, method="cluster")

    data['depth_cluster_plot'] = generate_clustering_plot(singlecell_df, tsne_results, cluster_labels, barcodes,
                                                          species_list, method="depth")

    if len(species_list) == 2:
        # Add an additional plot with cells colored by species for barnyard samples
        data['barnyard_cluster_plot'] = generate_clustering_plot(singlecell_df, tsne_results, cluster_labels, barcodes,
                                                                 species_list, method="species")

    helptext = [["Plots",
                 ["(left) Scatter plot of barcodes annotated as cells, colored by automatically computed clusters via graph clustering.",
                  "(right) Scatter plot of barcodes annotated as cells, colored by number of fragments in the barcode."]]]
    if is_barnyard:
        helptext[0][1] += ["(bottom) Scatter plot of barcodes annotated as cells, colored by species."]
    data['clustering_helptext'] = {'title': 'Cell Clustering', 'data': helptext}

    return data


def get_master_table(metadata, summary_data, species_list, is_barnyard, debug):
    data = {}
    if summary_data is None:
        return data
    metric_keys = [
        'annotated_cells',
        'cell_threshold',
        'median_fragments_per_cell',
        'median_fragments_per_noncell',
        'num_fragments',
        'total_usable_fragments',
        'frac_valid_barcode',
        'frac_valid_noncell',
        'frac_waste_overall_nondup',
        'waste_ratio_mito_cells',
        'r1_q30_bases_fract',
        'r2_q30_bases_fract',
        'bc_q30_bases_fract',
        'si_q30_bases_fract',
        'insert_twist_period',
        'insert_nucleosome_period',
        'frac_fragments_nfr',
        'frac_fragments_nuc',
        'frac_fragments_overlapping_targets',
        'frac_fragments_overlapping_tss',
        'frac_fragments_overlapping_dnase',
        'frac_fragments_overlapping_enhancer',
        'frac_fragments_overlapping_promoter',
        'frac_fragments_overlapping_blacklist',
        'frac_fragments_overlapping_peaks',
        'tss_enrichment_score',
        'ctcf_enrichment_score',
        'total_peaks_detected',
        'frac_cut_fragments_in_peaks',
        'bulk_total_library_complexity',
        'bulk_estimated_saturation',
        'frac_waste_duplicate',
        'frac_waste_mitochondrial',
        'median_per_cell_unique_fragments_at_{}_HQ_RPC'.format(RPC_10K),
        'median_per_cell_unique_fragments_at_{}_HQ_RPC'.format(RPC_30K),
        'median_per_cell_unique_fragments_at_{}_HQ_RPC'.format(RPC_50K),
        'median_per_cell_unique_fragments_at_{}_RRPC'.format(RPC_10K),
        'median_per_cell_unique_fragments_at_{}_RRPC'.format(RPC_30K),
        'median_per_cell_unique_fragments_at_{}_RRPC'.format(RPC_50K),
    ]
    if debug:
        metric_keys.append('estimated_gelbead_doublet_rate')
    if is_barnyard:
        metric_keys.extend([
            'observed_doublet_rate',
            'inferred_doublet_rate',
        ])
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list, debug=debug)
    rows = [[metric.category, metric.name, metric.value_string, metric.acceptable_string, metric.targeted_string]
            for metric in metrics]
    data['master_table'] = {
        "header": ["Category", "Metric", "Value", "Acceptable", "Targeted"],
        "rows": sorted(rows, key=lambda r: (r[0], r[1]))
    }

    data['summary_helptext'] = {'title': 'Summary Metrics', 'data': []}
    return data


def add_data(websummary_data, input_data):
    """Adds data to global dictionary, returns content"""
    if 'alarms' in input_data:
        websummary_data['alarms']['alarms'].extend(input_data['alarms'])
        del input_data['alarms']
    websummary_data.update(input_data)


def main(args, outs):
    reference = ReferenceManager(args.reference_path)
    species_list = reference.list_species()
    is_barnyard = len(species_list) > 1 and args.singlecell is not None

    singlecell_df = pd.read_csv(args.singlecell) if args.singlecell is not None else None

    summary_data = None
    if args.summary_results:
        with open(args.summary_results, 'r') as infile:
            summary_data = json.load(infile)

    # Pull up the correct template information
    template_path = os.path.dirname(os.path.abspath(__file__))
    template_file = os.path.join(template_path, '{}{}.html'.format(
        'barnyard' if is_barnyard else 'single',
        '_debug' if args.debug else ''))
    with open(template_file, 'r') as infile:
        template = infile.read()

    metadata = MetricAnnotations()
    websummary_data = {
        'alarms': {'alarms': []},
        'sample': {'id': args.sample_id,
                   'description': args.sample_desc,
                   'pipeline': "Cell Ranger ATAC"}
    }

    # Pull out all the general-purpose information
    add_data(websummary_data, get_hero_metric_data(metadata, summary_data, species_list, args.debug))
    add_data(websummary_data, get_pipeline_info(args, reference, args.debug))
    add_data(websummary_data, get_sequencing_info(metadata, summary_data, species_list, args.debug))
    add_data(websummary_data, get_cell_metrics_data(metadata, summary_data, species_list, singlecell_df,
                                                    args.excluded_barcodes, args.debug))
    add_data(websummary_data, get_clustering_plots(args.analysis, args.filtered_peak_bc_matrix, species_list,
                                                   singlecell_df, is_barnyard))
    add_data(websummary_data, get_insertsize_data(metadata, summary_data, singlecell_df, args.insert_sizes,
                                                  species_list, args.debug))
    add_data(websummary_data, get_targeting_data(metadata, summary_data, species_list, singlecell_df,
                                                 args.tss_relpos, args.ctcf_relpos, args.debug, DOWNSAMPLE_TARGETING))
    add_data(websummary_data, get_complexity_data(metadata, summary_data, args.bulk_complexity,
                                                  args.singlecell_complexity, species_list, args.debug))

    # For barnyard samples only
    if is_barnyard:
        add_data(websummary_data, get_barnyard_data(metadata, summary_data, species_list, singlecell_df,
                                                    args.debug, DOWNSAMPLE_BARNYARD))

    # For PD runs only
    if args.debug:
        add_data(websummary_data, get_peakcalling_data(metadata, summary_data, species_list, args.debug))
        add_data(websummary_data, get_wasted_data(metadata, summary_data, singlecell_df, species_list, args.debug))
        add_data(websummary_data, get_master_table(metadata, summary_data, species_list, is_barnyard, args.debug))

    # Modify the titles of plots to add consistent plot styling sample ID/descriptions
    for key, subdata in websummary_data.iteritems():
        if "layout" in subdata:
            subdata["layout"]["title"] += '<br><sup>Sample {} - {}</sup>'.format(args.sample_id, args.sample_desc)
            subdata["layout"]["hovermode"] = "closest"
            subdata["config"] = PLOT_CONFIG_KWARGS

    with open(outs.web_summary, 'w') as outfile:
        summarize.generate_html_summary(websummary_data, template, template_path, outfile)
