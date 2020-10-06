#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
# Runs alignment code
#
import itertools
import tenkit.align
import tenkit.fasta as tk_fasta
import tenkit.bam as tk_bam
import tenkit.reference
import martian
import math
from tools import ReferenceManager

__MRO__ = """
stage ALIGN_READS(
    in  map[] chunks,
    in  string aligner,
    in  string aligner_method,
    in  string reference_path,
    in  string read_group_sample,
    in  int    num_threads,
    out bam,
    src py     "stages/processing/align_reads",
) split using (
    in map chunk,
)
"""

# set a min read length -- BWA-MEM won't work below this
MIN_READ_LENGTH = 25

def split(args):
    """We just align each chunk independently -- joining will happen in the join step of SORT_READS"""

    # Pull some reads from fastq files -- bail out if it's less than 25bp
    fastq_tests = [x['read1'] for x in args.chunks]

    for fastq_test in fastq_tests:
        with open(fastq_test) as in_file:
            reader = tk_fasta.read_generator_fastq(in_file)
            for name, read, qual in itertools.islice(reader, 10):
                if len(read) < MIN_READ_LENGTH:
                    martian.alarm("BWA-MEM can't handle reads <25bp -- reads will be unmapped.")
                    continue

    # estimated amount of memory needed to process genome is 2x(num gigabases)+4GB
    ctg_mgr = ReferenceManager(args.reference_path)
    base_mem_in_gb = int(math.ceil(2 * ctg_mgr.get_vmem_est()))

    mem_in_gb = base_mem_in_gb + 4
    chunks = [{'chunk': x, '__threads': args.num_threads, '__mem_gb': mem_in_gb} for x in args.chunks]
    return {'chunks': chunks}

def main(args, outs):
    chunk = args.chunk

    if not chunk['reads_interleaved'] and (chunk['read1'] is None or chunk['read2'] is None):
        martian.throw("must supply a read1 and read2 when reads_interleave == False")

    if chunk['reads_interleaved']:
        reads = chunk['read1']
    else:
        reads = [chunk['read1']]
        if chunk['read2'] is not None:
            reads.append(chunk['read2'])

    a = tenkit.align.Aligner(reads, outs.default)
    aligner = args.aligner

    ref_fasta = tenkit.reference.get_fasta(args.reference_path)
    rg_string = chunk['read_group']
    read_group_header = tk_bam.make_rg_header(rg_string)
    a.output_alignment(aligner=aligner, aligner_params={'ref_fasta': ref_fasta, 'algorithm': args.aligner_method}, num_threads=martian.get_threads_allocation(), read_group_header=read_group_header)

def join(args, outs, chunk_defs, chunk_outs):
    outs.default = [chunk.default for chunk in chunk_outs]
