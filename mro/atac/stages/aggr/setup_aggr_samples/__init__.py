
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Determine normalization info per library
"""

# Imports
from __future__ import division
import pandas as pd
import numpy as np
from collections import Counter
import pickle
import json

from analysis.peaks import estimate_final_threshold
from tools import ReferenceManager, open_fragment_file
from constants import WINDOW_SIZE, PEAK_ODDS_RATIO
import martian

# MRO docstring
__MRO__ = """
stage SETUP_AGGR_SAMPLES(
    in  csv        aggr_csv,
    in  string     normalization,
    in  string     reference_path,
    out pickle     library_info,
    out map        gem_group_index,
    out json       gem_group_index_json,
    src py         "stages/aggr/setup_aggr_samples",
) split (
    in  int        n,
) using (
    volatile = strict,
)
"""

def split(args):
    """split into a chunk for each library in aggr csv, and define a unique gem group"""
    aggr_df = pd.read_csv(args.aggr_csv, sep=',')
    nchunks = len(aggr_df)

    ctg_mgr = ReferenceManager(args.reference_path)
    max_contig_len = max(ctg_mgr.get_contig_lengths().values())
    BYTES_PER_INT32_WITH_SAFETY = 5
    mem_gb = 2 * int(np.ceil(BYTES_PER_INT32_WITH_SAFETY * max_contig_len / 1024 / 1024 / 1024))

    return {'chunks': [{'n': group, '__mem_gb': mem_gb, '__vmem_gb': mem_gb + 6} for group in range(nchunks)], 'join': {'__mem_gb': 12}}

def main(args, outs):
    """Compute the depth and signal per library"""
    # read
    lib_id = args.n + 1
    aggr_df = pd.read_csv(args.aggr_csv, sep=',')
    library_info = {lib_id: {}}
    for label in aggr_df.columns.values.tolist():
        library_info[lib_id][label] = str(aggr_df.iloc[args.n][label])

    # if no normalization, don't waste compute
    if args.normalization is None:
        with open(outs.library_info, 'w') as f:
            pickle.dump(library_info, f)
        return

    # set ref properties
    ctg_mgr = ReferenceManager(args.reference_path)
    contig_lens = ctg_mgr.get_contig_lengths()
    max_contig_len = max(contig_lens.values())
    curr_chrom = None
    count_dict = Counter()
    chrom_len = 1
    half_window = WINDOW_SIZE // 2

    # traverse fragments file and count stats
    fragments_f = aggr_df.iloc[args.n]['fragments']
    Cuts = None
    special_normalization = (args.normalization in ["signal_mean", "signal_noise_threshold"])
    if special_normalization:
        Cuts = np.zeros(max_contig_len, dtype='int32')
    for chrom, start, stop, bc, dups in open_fragment_file(filename=fragments_f):
        if chrom != curr_chrom:
            curr_chrom = chrom
            if chrom not in contig_lens:
                martian.exit("fragment {}:{}-{} in {} is mapped to a contig not in the reference".format(chrom, start, stop, fragments_f))
            if special_normalization:
                count_dict += Counter(Cuts[i] for i in xrange(chrom_len) if Cuts[i] > 0)  # only traverse chrom len
                Cuts[:] = 0  # reset and reuse
                chrom_len = contig_lens[chrom]
        if special_normalization:
            Cuts[max(0, start - half_window): min(start + half_window + 1, chrom_len)] += 1
            Cuts[max(0, stop - half_window): min(stop + half_window + 1, chrom_len)] += 1
    if special_normalization:
        count_dict += Counter(Cuts[i] for i in xrange(chrom_len) if Cuts[i] > 0)  # only traverse chrom len

    scdf = pd.read_csv(library_info[lib_id]['cells'], sep=',')
    cell_mask = np.full(len(scdf), False)
    for species in ctg_mgr.list_species():
        cell_mask |= scdf['is_{}_cell_barcode'.format(species)] == 1
    library_info[lib_id]['total_fragments_per_cell'] = np.median(scdf[cell_mask]['total'] if 'total' in scdf[cell_mask].columns else scdf[cell_mask]['passed_filters'] + scdf[cell_mask]['duplicate'])
    library_info[lib_id]['unique_fragments_per_cell'] = np.median(scdf[cell_mask]['passed_filters'])

    # do peak calling fit on the count dict and get signal fit
    if args.normalization in ["signal_mean", "signal_noise_threshold"]:
        threshold, params = estimate_final_threshold(count_dict, PEAK_ODDS_RATIO)
        library_info[lib_id]['original_threshold'] = threshold
        library_info[lib_id]['signal_mean'] = 1 / params.p_signal

    # dump library info
    with open(outs.library_info, 'w') as f:
        pickle.dump(library_info, f)

def join(args, outs, chunk_defs, chunk_outs):
    """Compute per library subsample rate and the kind of normalization to do"""
    library_info = {}
    for chunk in chunk_outs:
        if chunk.library_info is not None:
            with open(chunk.library_info, 'r') as f:
                library_info.update(pickle.load(f))

    # trivial norm
    if args.normalization is None:
        for n in library_info:
            library_info[n]['kind'] = "unique"
            library_info[n]['rate'] = 1.0
        with open(outs.library_info, 'w') as f:
            pickle.dump(library_info, f)
    else:
        # other
        minval = 1e20  # arbitrarily large
        key = ""
        kind = ""
        if args.normalization == "unique":
            key = "unique_fragments_per_cell"
            kind = "unique"
        if args.normalization == "total":
            key = "total_fragments_per_cell"
            kind = "total"
        if args.normalization == "signal_noise_threshold":
            key = "original_threshold"
            kind = "unique"
        if args.normalization == "signal_mean":
            key = "signal_mean"
            kind = "unique"
        for n in library_info:
            library_info[n]['kind'] = kind
            minval = min(library_info[n][key], minval)
        for n in library_info:
            library_info[n]['rate'] = minval / library_info[n][key]
        # dump
        with open(outs.library_info, 'w') as f:
            pickle.dump(library_info, f)

    # derive this from library_info
    outs.gem_group_index = {"gem_group_index": {}}
    # populate the gem group index for legacy reasons where one could do an aggr between a lena metasample and another sample
    # in which case you'd need to remap the old groups in them metasample to new group ids. In ATAC, we are not allowing that behavior
    for lib_idx in library_info:
        outs.gem_group_index["gem_group_index"][lib_idx] = (library_info[lib_idx]['library_id'], 1)
    with open(outs.gem_group_index_json, 'w') as f:
        json.dump(outs.gem_group_index, f)
