"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Decide if cell calling should be done
"""

# Imports
import json

# MRO docstring
__MRO__ = """
stage DECIDE_CELL_CALLER(
    in  map  force_cells,
    in  csv  input_cell_barcodes,
    out bool disable_cell_calling,
    out json cell_calling_metrics,
    src py   "stages/reanalysis/decide_cell_caller",
)
"""

def main(args, outs):
    if args.force_cells is not None and args.input_cell_barcodes is not None:
        raise IOError("can't have force cells along with input choice of cell barcodes")

    # disable cell calling with input barcodes
    if args.input_cell_barcodes is not None:
        outs.disable_cell_calling = True

        # write out cell calling metrics
        counts = {}
        metrics = {}
        with open(args.input_cell_barcodes, 'r') as f:
            for line in f:
                items = line.strip("\n").rstrip(",").split(",")
                counts[items[0]] = len(items) - 1  # count bc per species
        for species in counts.keys():
            metrics["annotated_cells{}".format(species if len(counts) > 1 else "")] = counts[species]
        with open(outs.cell_calling_metrics, 'w') as f:
            json.dump(metrics, f, indent=4)
    else:
        # enable cell calling, whether force-cells is provided or not, unless input barcodes are provided
        outs.disable_cell_calling = False
        outs.cell_calling_metrics = None
