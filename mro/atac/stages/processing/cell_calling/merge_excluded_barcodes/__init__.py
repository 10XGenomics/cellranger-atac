"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Merges multiple JSON files with barcodes to be excluded from cell calling for
various reasons.
"""
from __future__ import absolute_import, division, print_function

import os
import json
from tenkit.safe_json import json_sanitize


__MRO__ = """
stage MERGE_EXCLUDED_BARCODES(
    in  json[] barcode_exclusions,
    out json   excluded_barcodes,
    src py     "stages/processing/cell_calling/merge_excluded_barcodes",
)
"""


def main(args, outs):
    exclusions = {}
    for filename in args.barcode_exclusions:
        if filename is None or not os.path.isfile(filename):
            continue
        with open(filename, "r") as infile:
            data = json.load(infile)
        reason = data["label"]
        for species, barcode_data in data["data"].iteritems():
            if species not in exclusions:
                exclusions[species] = {}
            for barcode, metric in barcode_data.iteritems():
                if barcode in exclusions[species]:
                    # This barcode was already excluded by another file
                    continue
                exclusions[species][barcode] = [reason, metric]

    with open(outs.excluded_barcodes, "w") as outfile:
        json.dump(json_sanitize(exclusions), outfile, indent=4, sort_keys=True)
