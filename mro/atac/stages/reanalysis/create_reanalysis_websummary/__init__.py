"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Produces the final customer-facing web_summary.
"""

from __future__ import division

from websummary import summarize

import json
import martian
import numpy as np
import h5py
import pandas as pd
import os

from plots.plotly_tools import PLOT_CONFIG_KWARGS
from plots.atac import generate_clustering_plot
from tools import ReferenceManager
from metrics import MetricAnnotations

import cellranger.matrix as cr_matrix
import cellranger.analysis.constants as cr_analysis_constants
import cellranger.analysis.graphclust as cr_graphclust

__MRO__ = """
stage CREATE_REANALYSIS_WEBSUMMARY(
    in  string reference_path,
    in  string barcode_whitelist,
    in  string sample_id,
    in  string sample_desc,
    in  csv    singlecell,
    in  json   summary_results,
    in  bool   debug,
    in  h5     filtered_peak_bc_matrix,
    in  h5     analysis,
    out html   web_summary,
    src py     "stages/reanalysis/create_reanalysis_websummary",
) using (
    mem_gb = 16,
)
"""

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
            data[key] = metric.gen_metric_dict(default_threshold='pass')

    return data

def get_pipeline_info(args, reference, debug):
    """Generates a table of general pipeline information.
    """
    data = {}
    metadata = reference.metadata

    rows = [
        ['Sample ID', args.sample_id],
        ['Sample description', args.sample_desc],
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
    return data

def get_clustering_plots(metadata, summary_data, analysis, filtered_peak_bc_matrix, species_list, singlecell_df, is_barnyard):
    data = {}
    if filtered_peak_bc_matrix is None or analysis is None:
        return data

    tsne_results = np.array(h5py.File(analysis, 'r')[cr_analysis_constants.ANALYSIS_H5_TSNE_GROUP]
                            ['_2']['transformed_tsne_matrix'])
    # N.B.: these cluster labels are 1-indexed
    cluster_labels = cr_graphclust.load_graphclust_from_h5(analysis).clusters
    barcodes = cr_matrix.CountMatrix.load_h5_file(filtered_peak_bc_matrix).bcs

    metric_keys = {
        'num_analysis_bcs',
    }

    # NOTE: HACK to register analysis barcodes without changing other stage code
    summary_data['num_analysis_bcs'] = len(barcodes)
    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list=species_list)
    data['analysis_table'] = {"rows": [[metric.name, metric.value_string] for metric in metrics]}
    helptext_metrics = metadata.gen_metric_helptext(metric_keys)

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
    data['clustering_helptext'] = {'title': 'Cell Clustering', 'data': helptext_metrics + helptext}

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
                   'pipeline': "Cell Ranger ATAC Renalyzer"}
    }

    singlecell_df = pd.read_csv(args.singlecell) if args.singlecell is not None else None

    add_data(websummary_data, get_hero_metric_data(metadata, summary_data, species_list, args.debug))

    add_data(websummary_data, get_pipeline_info(args, reference, args.debug))

    add_data(websummary_data, get_clustering_plots(metadata, summary_data, args.analysis,
             args.filtered_peak_bc_matrix, species_list, singlecell_df, is_barnyard))

    # Modify the titles of plots to add consistent plot styling sample ID/descriptions
    for key, subdata in websummary_data.iteritems():
        if "layout" in subdata:
            subdata["layout"]["title"] += '<br><sup>Sample {} - {}</sup>'.format(args.sample_id, args.sample_desc)
            subdata["layout"]["hovermode"] = "closest"
            subdata["config"] = PLOT_CONFIG_KWARGS

    with open(outs.web_summary, 'w') as outfile:
        summarize.generate_html_summary(websummary_data, template, template_path, outfile)
