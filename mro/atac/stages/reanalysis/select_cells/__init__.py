"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Select what file of cell barcodes should be used downstream in analysis
"""

# Imports
import cellranger.io as cr_io

# MRO docstring
__MRO__ = """
stage SELECT_CELLS(
    in  bool selector,
    in  csv  cell_barcodes1,
    in  csv  cell_barcodes2,
    in  json metrics1,
    in  json metrics2,
    out csv  cell_barcodes,
    out json metrics,
    src py   "stages/reanalysis/select_cells",
)
"""

def main(args, outs):
    if args.selector:
        if args.cell_barcodes1 is None:
            raise IOError("Input barcodes file 1 is not present")
        cr_io.copy(args.cell_barcodes1, outs.cell_barcodes)
        cr_io.copy(args.metrics1, outs.metrics)
    else:
        if args.cell_barcodes2 is None:
            raise IOError("Input barcodes file 2 is not present")
        cr_io.copy(args.cell_barcodes2, outs.cell_barcodes)
        cr_io.copy(args.metrics2, outs.metrics)
