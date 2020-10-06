
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

copies inputs to outputs
"""

# Imports
import martian
import cellranger.io as cr_io

# MRO docstring
__MRO__ = """
stage COPY_INPUT_OUTPUT(
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  csv        aggr_csv,
    out tsv.gz     fragments,
    out tsv.gz.tbi fragments_index,
    out csv        aggr_csv,
    src py         "stages/common/copy_input_output",
) split (
)
"""

def split(args):
    return {'chunks': [], 'join': {}}

def main(args, outs):
    martian.throw("main not implemented")

def join(args, outs, chunk_defs, chunk_outs):
    for infile, outfile in zip([args.fragments, args.fragments_index, args.aggr_csv], [outs.fragments, outs.fragments_index, outs.aggr_csv]):
        if infile is None:
            outfile = infile
        else:
            cr_io.copy(infile, outfile)
