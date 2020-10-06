"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Collates basic information about reads from a BAM file
"""
from __future__ import division

import os

import martian

import tenkit.bio_io as tk_io
import tenkit.fasta as tk_fasta
import tenkit.safe_json
import tenkit.summary_manager as tk_summary
from tenkit.stats import robust_divide

from barcodes import get_read_raw_barcode, get_read_barcode_qual
from tools import create_bam_infile

__MRO__ = """
stage GENERATE_SEQUENCING_METRICS(
    in  bam[]  input,
    out txt    misc_sm,
    out json   summary,
    src py     "stages/metrics/generate_bulk_metrics",
) split using (
    in  bam    chunk_bam,
)
"""


def combine_summary_managers(sm_pickles):
    """Load in a list of pickled SummaryManager objects and merge them into a single object
    """
    summary_manager = None
    for filename in sm_pickles:
        tmp_sm = tk_summary.SummaryManager.load(filename)
        if summary_manager is None:
            summary_manager = tmp_sm
        else:
            summary_manager = summary_manager.combine(tmp_sm)
    return summary_manager


def split(args):
    chunk_defs = [{'chunk_bam': x, '__mem_gb': 3} for x in args.input if os.path.isfile(x)]
    return {'chunks': chunk_defs, 'join': {'__mem_gb': 8}}


def join(args, outs, chunk_defs, chunk_outs):
    """Joins the various chunk outputs and computes
    further summary metrics based on the merged metrics."""
    martian.log_info("Combining miscellaneous summary managers")
    misc_sm = combine_summary_managers([chunk.misc_sm for chunk in chunk_outs])

    martian.log_info("Computing summary metrics")
    compute_summary_metrics(misc_sm)

    with open(outs.summary, 'w') as outfile:
        outfile.write(tenkit.safe_json.safe_jsonify(misc_sm.get_summarizer('metrics').dict, pretty=True))


def main(args, outs):
    bam_in = create_bam_infile(args.chunk_bam)
    references = bam_in.references
    misc_sm = compute_basic_stats(bam_in, references)
    misc_sm.save(outs.misc_sm)


def compute_basic_stats(bam_in, references):
    """Gathers basic information about the reads in a single chunk.
    """

    # Initialize all the summarizers and load them into a manager
    metrics = tk_summary.DictionaryDistribution()
    misc_sm = tk_summary.SummaryManager()
    misc_sm.add_summarizer(metrics, 'metrics')

    for info in read_pair_info_iter(bam_in, references):
        metrics.add('num_reads', 2)

        r1_seq_len = info['r1_seq_len']
        r2_seq_len = info['r2_seq_len']

        r1_q30_bases = info['r1_q30_bases']
        if not (r1_q30_bases is None or r1_seq_len is None or r1_seq_len == 0):
            metrics.add('r1_q30_bases', r1_q30_bases)
            metrics.add('r1_tot_bases', r1_seq_len)

        r2_q30_bases = info['r2_q30_bases']
        if not (r2_q30_bases is None or r2_seq_len is None or r2_seq_len == 0):
            metrics.add('r2_q30_bases', r2_q30_bases)
            metrics.add('r2_tot_bases', r2_seq_len)

        si_q30_bases = info['si_q30_bases']
        si_tot_bases = info['sample_index_len']
        if not (si_q30_bases is None or si_tot_bases is None or si_tot_bases == 0):
            metrics.add('si_q30_bases', si_q30_bases)
            metrics.add('si_tot_bases', si_tot_bases)

        bc_q30_bases = info['bc_q30_bases']
        bc_tot_bases = info['10X_raw_bc_len']
        if not (bc_q30_bases is None or bc_tot_bases is None or bc_tot_bases == 0):
            metrics.add('bc_q30_bases', bc_q30_bases)
            metrics.add('bc_tot_bases', bc_tot_bases)

    return misc_sm


def compute_summary_metrics(misc_sm):
    """Called in the join step to extract summary metrics from the pooled summarizer objects"""
    metrics = misc_sm.get_summarizer('metrics')

    metrics['r1_q30_bases_fract'] = robust_divide(metrics['r1_q30_bases'], metrics['r1_tot_bases'])
    metrics['r2_q30_bases_fract'] = robust_divide(metrics['r2_q30_bases'], metrics['r2_tot_bases'])
    metrics['si_q30_bases_fract'] = robust_divide(metrics['si_q30_bases'], metrics['si_tot_bases'])
    metrics['bc_q30_bases_fract'] = robust_divide(metrics['bc_q30_bases'], metrics['bc_tot_bases'])

    return metrics


def read_pair_info_iter(bam_in, references):
    """Takes a read-sorted bam file along with references appropriate to that BAM file and reports various information
    about each read pair in the file.

    Does not support single-end sequencing data.
    """
    while True:
        # First read in the pair, skipping any secondary alignments
        try:
            read1 = bam_in.next()
            while read1.is_secondary:
                read1 = bam_in.next()
            read2 = bam_in.next()
            while read2.is_secondary:
                read2 = bam_in.next()
        except StopIteration:
            break

        name1 = read1.qname.split()[0]
        name2 = read2.qname.split()[0]
        assert (name1 == name2)

        # Swap reads so that read 1 is first
        if read1.reference_id == read2.reference_id and read2.reference_start < read1.reference_start:
            read1, read2 = read2, read1

        read_pair_info = {
            'name': name1,
            'r1_seq_len': 0 if read1.seq is None else len(read1.seq),
            'r2_seq_len': 0 if read2.seq is None else len(read2.seq),
            'r1_q30_bases': tk_fasta.get_bases_qual(read1.qual, 30),
            'r2_q30_bases': tk_fasta.get_bases_qual(read2.qual, 30),
        }

        read_raw_bc = get_read_raw_barcode(read1, '')
        read_bc_qual = get_read_barcode_qual(read1, '')

        # Note these don't have overridable defaults and so may be None
        read_sample_index = tk_io.get_read_sample_index(read1)
        read_sample_index = '' if read_sample_index is None else read_sample_index
        read_sample_index_qual = tk_io.get_read_sample_index_qual(read1)
        read_sample_index_qual = '' if read_sample_index_qual is None else read_sample_index_qual

        read_pair_info['10X_raw_bc_len'] = len(read_raw_bc)
        read_pair_info['sample_index_len'] = len(read_sample_index)
        read_pair_info['si_q30_bases'] = tk_fasta.get_bases_qual(read_sample_index_qual, 30)
        read_pair_info['bc_q30_bases'] = tk_fasta.get_bases_qual(read_bc_qual, 30)

        yield read_pair_info
