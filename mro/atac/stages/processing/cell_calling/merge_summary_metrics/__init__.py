"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Merges multiple JSON files with summary metrics produced from different stages into a
single output file.
"""
from __future__ import absolute_import, division, print_function

import os
import json
from tenkit.safe_json import json_sanitize


__MRO__ = """
stage MERGE_SUMMARY_METRICS(
    in  json[] summary_jsons,
    out json   merged_summary,
    src py     "stages/processing/cell_calling/merge_summary_metrics",
)
"""


def main(args, outs):
    summary = {}
    for filename in args.summary_jsons:
        if filename is None or not os.path.isfile(filename):
            continue
        with open(filename, "r") as infile:
            summary.update(json.load(infile))

    with open(outs.merged_summary, "w") as outfile:
        json.dump(json_sanitize(summary), outfile, indent=4, sort_keys=True)
