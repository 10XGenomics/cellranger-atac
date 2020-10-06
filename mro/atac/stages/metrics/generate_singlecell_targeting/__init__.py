"""
Generate targeting metrics for individual barcodes

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division
import numpy as np
import pickle
from tools import ReferenceManager
from tools.io import get_counts_by_barcode
from collections import Counter

__MRO__ = """
stage GENERATE_SINGLECELL_TARGETING(
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  bed        peaks,
    in  string     reference_path,
    out csv        singlecell,
    out json       summary,
    out csv        tss_relpos,
    out csv        ctcf_relpos,
    src py         "stages/metrics/generate_singlecell_targeting",
) split (
    in  string     contig,
    out int        read_count,
    out pickle     target_counts_by_barcode,
    out pickle     chunk_tss,
    out pickle     chunk_ctcf,
) using (
    mem_gb   = 6,
    volatile = strict,
)
"""


def split(args):
    if args.fragments is None:
        return {'chunks': [], 'join': {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    chunks = []
    for contig in all_contigs:
        chunks.append({'contig': contig,
                       '__mem_gb': 6})

    return {'chunks': chunks, 'join': {'__mem_gb': 6}}


def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    if args.fragments is None:
        outs.read_count = None
        outs.singlecell = None
        outs.tss_relpos = None
        outs.ctcf_relpos = None
        return

    read_count = 0
    target_counts_by_barcode = {}
    tss_relpos = Counter()
    ctcf_relpos = Counter()
    for chunk_in, chunk_out in zip(chunk_defs, chunk_outs):
        read_count += chunk_out.read_count

        with open(chunk_out.target_counts_by_barcode, 'r') as infile:
            chunk_counts = pickle.load(infile)
        for barcode, barcode_counts in chunk_counts.iteritems():
            if barcode not in target_counts_by_barcode:
                target_counts_by_barcode[barcode] = Counter()
            target_counts_by_barcode[barcode] += barcode_counts

        with open(chunk_out.chunk_tss, 'r') as infile:
            chunk_tss = pickle.load(infile)
        tss_relpos += chunk_tss

        with open(chunk_out.chunk_ctcf, 'r') as infile:
            chunk_ctcf = pickle.load(infile)
        ctcf_relpos += chunk_ctcf

    keys = ["{region}_fragments".format(region=reg)
            for reg in ["TSS", "DNase_sensitive_region", "enhancer_region",
                        "promoter_region", "on_target", "blacklist_region", "peak_region"]] + ["peak_region_cutsites"]

    with open(outs.singlecell, 'w') as outfile:
        outfile.write("barcode,")
        outfile.write(",".join(keys))
        outfile.write("\n")
        for barcode in sorted(target_counts_by_barcode.keys()):
            outfile.write("{},".format(barcode))
            outfile.write(",".join([str(target_counts_by_barcode[barcode][key]) for key in keys]))
            outfile.write("\n")

    # Normalize the TSS profile to the total # of TSS sites
    ref_manager = ReferenceManager(args.reference_path)
    if ref_manager.tss_track is None:
        outs.tss_relpos = None
    else:
        keys = range(-2000, 2001)
        with open(ref_manager.tss_track, 'r') as infile:
            tss_count = len(infile.readlines())
        try:
            tss_norm_factor = 10e6 / (tss_count * read_count)
        except ZeroDivisionError:
            tss_norm_factor = np.nan

        with open(outs.tss_relpos, 'w') as outfile:
            outfile.write(",".join([str(k) for k in keys]))
            outfile.write("\n")
            outfile.write(",".join([str(tss_relpos[key] * tss_norm_factor) for key in keys]))
            outfile.write("\n")

    if ref_manager.ctcf_track is None:
        outs.ctcf_relpos = None
    else:
        keys = range(-250, 251)
        with open(outs.ctcf_relpos, 'w') as outfile:
            outfile.write(",".join([str(k) for k in keys]))
            outfile.write("\n")
            outfile.write(",".join([str(ctcf_relpos[key]) for key in keys]))
            outfile.write("\n")


def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    if args.fragments is None:
        return

    read_count, target_counts, tss_relpos, ctcf_relpos = get_counts_by_barcode(args.reference_path, args.peaks,
                                                                               args.fragments,
                                                                               fragments_index=args.fragments_index,
                                                                               contig=args.contig)
    outs.read_count = read_count
    with open(outs.target_counts_by_barcode, 'w') as outfile:
        pickle.dump(target_counts, outfile)
    with open(outs.chunk_tss, 'w') as outfile:
        pickle.dump(tss_relpos, outfile)
    with open(outs.chunk_ctcf, 'w') as outfile:
        pickle.dump(ctcf_relpos, outfile)
