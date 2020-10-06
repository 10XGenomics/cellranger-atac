"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Sorts the input BAM file by position, using the 5' position tag of each read as the primary key and the
alignment position as a secondary key.

Using the 5' position tag improves duplicate marking by "folding down" softclipped bases at the 5' ends of reads.
"""

import os.path

import martian

from constants import SELF_FIVE_PRIME_POS_TAG
from tools.io import hierarchical_merge_bam, sort_bam

__MRO__ = """
stage SORT_READS_BY_POS(
    in  bam[]  input,
    out bam    tagsorted_bam,
    src py     "stages/processing/sort_reads_by_pos",
) split (
    in  bam chunk_input,
)
"""


def split(args):
    chunk_defs = [{"chunk_input": x} for x in args.input if os.path.isfile(x)]
    return {"chunks": chunk_defs, "join": {"__threads": 4}}


def main(args, outs):
    """Use samtools to sort each chunk using the 5-prime position tag, using (by default) the alignment position
    as a secondary key.
    """
    args.coerce_strings()
    bam_prefix, _ = os.path.splitext(outs.tagsorted_bam)
    sort_bam(args.chunk_input, "{}.bam".format(bam_prefix), SELF_FIVE_PRIME_POS_TAG,
             martian.get_threads_allocation())


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.tagsorted_bam) for chunk in chunk_outs]
    hierarchical_merge_bam(input_bams, outs.tagsorted_bam, SELF_FIVE_PRIME_POS_TAG,
                           martian.get_threads_allocation())
