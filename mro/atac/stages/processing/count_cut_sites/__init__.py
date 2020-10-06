"""
Detect open chromatin regions in the genome as measured from ATAC-seq fragment ends.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
from __future__ import division

import pickle
from tools import ReferenceManager, combine_csv, parsed_fragments_from_contig
from collections import Counter
import numpy as np
from constants import WINDOW_SIZE

__MRO__ = """
stage COUNT_CUT_SITES(
    in  path       reference_path,
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    out bedgraph   cut_sites,
    out pickle     count_dict,
    src py         "stages/processing/count_cut_sites",
) split (
    in  string     contig,
)
"""

def write_chrom_bedgraph(chrom, contig_len, Cuts, out_bedgraph, min_counts=1):
    last_count = None
    last_start = None
    last_pos = None

    with open(out_bedgraph, 'w') as outfile_bg:
        for pos in xrange(contig_len):
            # NB that the output dict always keeps track of windowed counts
            count = Cuts[pos]
            if count == 0:
                continue
            if last_count == count:
                # We can continue a section of the bed graph
                last_pos = pos
            else:
                if last_count >= min_counts:
                    outfile_bg.write('{}\t{}\t{}\t{}\n'.format(chrom, last_start, last_pos + 1, last_count))
                last_start = pos
                last_pos = pos
                last_count = count
        if last_count >= min_counts:
            outfile_bg.write('{}\t{}\t{}\t{}\n'.format(chrom, last_start, last_pos + 1, last_count))

def split(args):
    if args.fragments is None:
        return {'chunks': [], 'join': {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)
    contig_len = ctg_mgr.get_contig_lengths()
    BYTES_PER_INT32_WITH_SAFETY = 5

    chunks = []
    for contig in all_contigs:
        chunks.append({'contig': contig,
                       '__mem_gb': int(np.ceil(BYTES_PER_INT32_WITH_SAFETY * contig_len[contig] / 1024 / 1024 / 1024))})

    return {'chunks': chunks, 'join': {'__mem_gb': 5}}

def join(args, outs, chunk_defs, chunk_outs):
    '''Merge counts and bedgraph'''
    if args.fragments is None:
        outs.count_dict = None
        outs.cut_sites = None
        return

    # merge bedgraphs in contig order in chunks
    chunk_bedgraphs = [chunk_out.cut_sites for chunk_out in chunk_outs if chunk_out.cut_sites is not None]
    combine_csv(chunk_bedgraphs, outs.cut_sites, header_lines=0)

    # Merge count dicts
    count_dict = Counter()
    for chunk_out in chunk_outs:
        with open(chunk_out.count_dict, 'r') as f:
            count_dict += pickle.load(f)
    with open(outs.count_dict, 'w') as f:
        pickle.dump(count_dict, f)

def main(args, outs):
    '''Find cut sites on a per chromosome basis and write out a bedgraph'''
    if args.fragments is None:
        outs.count_dict = None
        outs.cut_sites = None
        return

    ctg_mgr = ReferenceManager(args.reference_path)
    contig_len = ctg_mgr.get_contig_lengths()
    chrom_len = contig_len[args.contig]
    half_window = WINDOW_SIZE // 2
    Cuts = np.zeros(chrom_len, dtype='int32')

    # find windowed cut sites
    for _, start, stop, _, _ in parsed_fragments_from_contig(contig=args.contig, filename=args.fragments, index=args.fragments_index):
        Cuts[max(0, start - half_window): min(start + half_window + 1, chrom_len)] += 1
        Cuts[max(0, stop - half_window): min(stop + half_window + 1, chrom_len)] += 1

    # get count dict
    count_dict = Counter(v for v in Cuts if v > 0)
    with open(outs.count_dict, 'w') as count_dict_out:
        pickle.dump(count_dict, count_dict_out)

    # write bedgraph of * windowed cutsites *
    if len(count_dict):
        write_chrom_bedgraph(args.contig, chrom_len, Cuts, outs.cut_sites)
    else:
        outs.cut_sites = None
