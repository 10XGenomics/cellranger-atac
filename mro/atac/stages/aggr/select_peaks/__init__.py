"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Select what file of  peaks should be used downstream in analysis
"""

# Imports
import cellranger.io as cr_io

# MRO docstring
__MRO__ = """
stage SELECT_CELLS(
    in  bool selector,
    in  csv  peaks1,
    in  csv  peaks2,
    out csv  peaks,
    src py   "stages/reanalysis/select_peaks",
)
"""

def main(args, outs):
    if args.selector:
        if args.peaks1 is None:
            raise IOError("Input peaks file 1 is not present")
        cr_io.copy(args.peaks1, outs.peaks)
    else:
        if args.peaks2 is None:
            raise IOError("Input peaks file 2 is not present")
        cr_io.copy(args.peaks2, outs.peaks)
