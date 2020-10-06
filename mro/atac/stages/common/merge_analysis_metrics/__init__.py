
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

aggregates metrics
"""

# Imports
import json
import pickle

import martian
import numpy as np
import tenkit.safe_json
from metrics import MetricAnnotations
from tools.io import write_dict_to_csv
from tools import ReferenceManager

# MRO docstring
__MRO__ = """
stage MERGE_ANALYSIS_METRICS(
    in  string reference_path,
    in  pickle library_info,
    in  json[] metrics,
    in  json   normalization_metrics,
    in  json   peak_calling_metrics,
    in  json   aggregation_metrics,
    out json   metrics,
    out csv    metrics_csv,
    src py     "stages/common/merge_analysis_metrics",
)
"""

def main(args, outs):
    metrics = {}
    for fname in args.metrics:
        if fname is not None:
            with open(fname, 'r') as f:
                metrics.update(json.load(f))

    # Normalize "NaN" values
    for key in metrics:
        value = metrics[key]
        if str(value) == 'NaN' or (isinstance(value, float) and np.isnan(value)):
            metrics[key] = None

    # add version info
    metrics['cellranger-atac_version'] = martian.get_pipelines_version()

    if len(metrics) > 0:
        martian.log_info('Writing out summary_metrics')
        with open(outs.metrics, 'w') as outfile:
            outfile.write(tenkit.safe_json.safe_jsonify(metrics, pretty=True))

    # compile summary.csv metrics
    # load library info and fake libraries as species
    metric_registry = MetricAnnotations()
    metrics_csv_dict = {}
    if args.library_info is not None:
        with open(args.library_info, 'r') as f:
            library_info = pickle.load(f)
        library_list = [library_info[n]['library_id'] for n in library_info.keys()]
        metrics_csv_dict.update(metric_registry.compile_summary_metrics(metrics, species_list=library_list))

    # load species level metrics
    ctg_mgr = ReferenceManager(args.reference_path)
    metrics_csv_dict.update(metric_registry.compile_summary_metrics(metrics, species_list=ctg_mgr.list_species()))
    write_dict_to_csv(outs.metrics_csv, metrics_csv_dict, sort=True)
