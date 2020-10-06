"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Merges data from other stages into comprehensive outputs
"""

from __future__ import division

import json
import numpy as np

import martian

import tenkit
import tenkit.hdf5
import tenkit.safe_json
from tools.io import write_dict_to_csv
from tools import ReferenceManager

from metrics import MetricAnnotations

__MRO__ = """
stage SUMMARIZE_REPORTS_SINGLECELL(
    in  string reference_path,
    in  json   complexity_summary,
    in  json   cell_calling_summary,
    in  json   peak_results,
    in  json   basic_results,
    in  json   error_results_summary,
    in  json   insert_summary,
    in  json   singlecell_results,
    in  json   contam_results,
    in  json   downsample_info,
    in  json   enrichment_results,
    out json   analysis_params,
    out json   summary,
    out csv    summary_csv,
    src py     "stages/reporter/summarize_reports_singlecell",
) using (
    mem_gb = 4,
)
"""


def split(args):
    raise NotImplementedError('No split')


def join(args, outs, chunk_defs, chunk_outs):
    raise NotImplementedError('No split')


def write_analysis_parameters(analysis_params_outfn):
    with open(analysis_params_outfn, 'w') as analysis_params_out:
        analysis_params = {
            'analysis_version': martian.get_pipelines_version(),
            # Dropping meowmix version -- we're moving to putting special reference datasets into main repo
            'meowmix_version': "99.9.9",
            # Lena needs this set, even though we're not trimming
            'lead_trim': 0,
        }
        analysis_params_out.write(tenkit.safe_json.safe_jsonify(analysis_params))


def simple_load_metrics(summary_metrics, metrics_fn):
    with open(metrics_fn, 'r') as infile:
        metrics = json.load(infile)
    summary_metrics.update(metrics)
    summary_metrics['cellranger-atac_version'] = martian.get_pipelines_version()
    return summary_metrics


def main(args, outs):
    reference = ReferenceManager(args.reference_path)

    martian.log_info('Writing analysis parameters')
    write_analysis_parameters(outs.analysis_params)

    martian.log_info('Initializing summary metrics')
    summary_metrics = {}
    summary_metrics = simple_load_metrics(summary_metrics, args.basic_results)

    if args.singlecell_results is not None:
        martian.log_info('Loading single cell results')
        summary_metrics = simple_load_metrics(summary_metrics, args.singlecell_results)

    if args.insert_summary is not None:
        martian.log_info('Loading insert summary')
        summary_metrics = simple_load_metrics(summary_metrics, args.insert_summary)

    if args.complexity_summary is not None:
        martian.log_info('Loading complexity summary')
        summary_metrics = simple_load_metrics(summary_metrics, args.complexity_summary)

    if args.error_results_summary is not None:
        martian.log_info('Loading error summary')
        summary_metrics = simple_load_metrics(summary_metrics, args.error_results_summary)

    if args.downsample_info is not None:
        martian.log_info('Loading downsampling information')
        summary_metrics = simple_load_metrics(summary_metrics, args.downsample_info)

    if args.contam_results is not None:
        martian.log_info('Loading contamination results')
        summary_metrics = simple_load_metrics(summary_metrics, args.contam_results)

    if args.peak_results is not None:
        martian.log_info('Loading peak results')
        summary_metrics = simple_load_metrics(summary_metrics, args.peak_results)

    if args.enrichment_results is not None:
        martian.log_info('Loading TSS and CTCF scores')
        summary_metrics = simple_load_metrics(summary_metrics, args.enrichment_results)

    if args.cell_calling_summary is not None:
        martian.log_info('Loading cell calling parameters')
        summary_metrics = simple_load_metrics(summary_metrics, args.cell_calling_summary)

    # Normalize "NaN" values
    for key in summary_metrics:
        value = summary_metrics[key]
        if str(value) == 'NaN' or (isinstance(value, float) and np.isnan(value)):
            summary_metrics[key] = None

    if reference.metadata:
        # If we have reference metadata - copy over the data to summary.json
        for (key, value) in reference.metadata.items():
            summary_metrics["reference_" + key] = value

    martian.log_info('Writing out summary_metrics')
    with open(outs.summary, 'w') as outfile:
        outfile.write(tenkit.safe_json.safe_jsonify(summary_metrics, pretty=True))

    # compile summary.csv metrics
    metric_registry = MetricAnnotations()

    species_list = reference.list_species()
    summary_csv_dict = metric_registry.compile_summary_metrics(summary_metrics, species_list=species_list)
    write_dict_to_csv(outs.summary_csv, summary_csv_dict, sort=True)
