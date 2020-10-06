"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Process a pipestance outs directory and harvest secondary analysis outputs for validation
"""

import os

__MRO__ = """
stage HARVEST_PIPESTANCE(
    in  path   pipestance,
    out bam    possorted_bam,
    out path   analysis,
    out tsv.gz fragments,
    src py   "stages/validate/harvest_pipestance",
)
"""

def main(args, outs):
    """Run the main code. Note that the outs structure defines the processing below"""

    if not os.path.exists(args.pipestance):
        raise IOError("pipestance {} does not exist".format(args.pipestance))

    pipe_outs = os.path.join(args.pipestance, "outs")
    if not os.path.exists(pipe_outs):
        raise IOError("outs directory not found")

    possorted_bam = os.path.join(pipe_outs, "possorted_bam.bam")
    possorted_bam_index = os.path.join(pipe_outs, "possorted_bam.bam.bai")
    if not os.path.exists(possorted_bam) or not os.path.exists(possorted_bam_index):
        raise IOError("position sorted bam or its index does not exist in the pipestance outs")

    outs.possorted_bam = possorted_bam

    analysis = os.path.join(pipe_outs, "analysis")
    if not os.path.exists(analysis):
        raise IOError("no analysis files found")
    outs.analysis = analysis

    fragments = os.path.join(pipe_outs, "fragments.tsv.gz")
    if not os.path.exists(fragments):
        raise IOError("fragments not found")
    outs.fragments = fragments
