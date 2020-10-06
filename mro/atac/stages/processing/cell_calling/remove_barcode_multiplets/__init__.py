"""
Identifies likely barcode multiplets.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
from __future__ import absolute_import, division, print_function

from barcodes import query_barcode_subsequences, split_barcode, merge_barcode
from tools import ReferenceManager, grouped_fragments_from_contig
import numpy as np
import json
from collections import Counter
import gzip

__MRO__ = """
stage REMOVE_BARCODE_MULTIPLETS(
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  string     reference_path,
    in  string     barcode_whitelist,
    in  json       barcode_counts,
    out json       barcode_multiplets,
    out json       barcode_multiplets_summary,
    src py         "stages/processing/cell_calling/remove_barcode_multiplets",
) split (
    in  string     gem_group,
    in  string     contig,
    out npy.gz     part_a_linkage_matrix,
    out npy.gz     part_b_linkage_matrix,
) using (
    mem_gb = 4,
    volatile = strict,
)
"""

MAXIMUM_FRAGMENT_SIZE = 800
MAXIMUM_POSITION_SIZE = 2500
MINIMUM_COUNT = 100
MAXIMUM_PILEUP = 20

SELF_SIGNAL_THRESHOLD_MULTIPLIER = 0.4


def split(args):
    if args.fragments is None:
        return {"chunks": [], "join": {}}

    with open(args.barcode_counts, "r") as infile:
        barcode_counts = Counter(json.load(infile))

    valid_barcodes = barcode_counts.keys()
    part_a_seqs, part_c_seqs, part_b_seqs, gem_group_seqs = query_barcode_subsequences(
        valid_barcodes
    )

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    chunks = []
    for gem_group in gem_group_seqs:
        for contig in all_contigs:
            chunks.append(
                {
                    "contig": contig,
                    "gem_group": gem_group,
                    "__mem_gb": 4,
                }
            )

    return {"chunks": chunks, "join": {"__mem_gb": 16}}


def join(args, outs, chunk_defs, chunk_outs):
    if args.fragments is None:
        outs.barcode_multiplets = None
        outs.barcode_multiplets_summary = None
        return

    with open(args.barcode_counts, "r") as infile:
        barcode_counts = Counter(json.load(infile))

    valid_barcodes = barcode_counts.keys()
    part_a_seqs, part_c_seqs, part_b_seqs, gem_group_seqs = query_barcode_subsequences(
        valid_barcodes
    )

    part_a_count = max([len(part_a_seqs[c]) for c in part_c_seqs])
    part_b_count = max([len(part_b_seqs[c]) for c in part_c_seqs])
    part_c_count = len(part_c_seqs)

    index_by_part_a = {
        part_c: {
            part_a: i for i, part_a in enumerate(part_a_seqs[part_c])
        } for part_c in part_c_seqs
    }
    index_by_part_b = {
        part_c: {
            part_b: i for i, part_b in enumerate(part_b_seqs[part_c])
        } for part_c in part_c_seqs
    }
    index_by_part_c = {
        part_c: i for i, part_c in enumerate(part_c_seqs)
    }

    part_a_linkage_matrix = np.zeros((part_c_count, part_b_count, part_a_count, part_a_count), dtype=np.uint32)
    part_b_linkage_matrix = np.zeros((part_c_count, part_a_count, part_b_count, part_b_count), dtype=np.uint32)

    # Search for contaminants as barcodes with higher similarity to a major barcode
    # with some mimimum signal than self-similarity.
    barcode_multiplets = {}

    # group chunks by gem group and aggregate across contigs for post-processing
    for gem_group_seq in gem_group_seqs:
        part_a_linkage_matrix[:, :, :, :] = 0
        part_b_linkage_matrix[:, :, :, :] = 0

        for chunk_in, chunk_out in zip(chunk_defs, chunk_outs):
            if gem_group_seq != chunk_in.gem_group:
                continue

            # aggregate across contigs
            infile = gzip.GzipFile(chunk_out.part_a_linkage_matrix, 'r')
            part_a_linkage_matrix += np.load(infile)
            infile.close()

            infile = gzip.GzipFile(chunk_out.part_b_linkage_matrix, 'r')
            part_b_linkage_matrix += np.load(infile)
            infile.close()

        for major_barcode, count in barcode_counts.iteritems():
            if count < MINIMUM_COUNT:
                continue
            part_a, part_c, part_b, gem_group = split_barcode(major_barcode, return_gg=True)
            if gem_group != gem_group_seq:
                continue

            part_a_index = index_by_part_a[part_c][part_a]
            part_b_index = index_by_part_b[part_c][part_b]
            part_c_index = index_by_part_c[part_c]

            for other_part_a in part_a_seqs[part_c]:
                if other_part_a == part_a:
                    continue
                minor_barcode = merge_barcode(other_part_a, part_c, part_b, gem_group)

                other_part_a_index = index_by_part_a[part_c][other_part_a]
                self_signal = part_a_linkage_matrix[part_c_index, part_b_index, other_part_a_index, other_part_a_index]
                major_signal = part_a_linkage_matrix[part_c_index, part_b_index, other_part_a_index, part_a_index]
                if major_signal > (self_signal * SELF_SIGNAL_THRESHOLD_MULTIPLIER):
                    if minor_barcode not in barcode_multiplets:
                        barcode_multiplets[minor_barcode] = major_barcode
                    else:
                        old_major = barcode_multiplets[minor_barcode]
                        old_a, _, _ = split_barcode(old_major)
                        old_a_index = index_by_part_a[part_c][old_a]
                        old_signal = part_a_linkage_matrix[part_c_index, part_b_index, other_part_a_index, old_a_index]
                        if major_signal > old_signal:
                            barcode_multiplets[minor_barcode] = major_barcode

            for other_part_b in part_b_seqs[part_c]:
                if other_part_b == part_b:
                    continue
                minor_barcode = merge_barcode(part_a, part_c, other_part_b, gem_group)

                other_part_b_index = index_by_part_b[part_c][other_part_b]
                self_signal = part_b_linkage_matrix[part_c_index, part_a_index, other_part_b_index, other_part_b_index]
                major_signal = part_b_linkage_matrix[part_c_index, part_a_index, other_part_b_index, part_b_index]
                if major_signal > (self_signal * SELF_SIGNAL_THRESHOLD_MULTIPLIER):
                    if minor_barcode not in barcode_multiplets:
                        barcode_multiplets[minor_barcode] = major_barcode
                    else:
                        old_major = barcode_multiplets[minor_barcode]
                        _, _, old_b = split_barcode(old_major)
                        old_b_index = index_by_part_b[part_c][old_b]
                        old_signal = part_b_linkage_matrix[
                            part_c_index, part_a_index, other_part_b_index, old_b_index]
                        if major_signal > old_signal:
                            barcode_multiplets[minor_barcode] = major_barcode

    # Post-screen the contaminants for pairs that are linked to each other.  In that
    # case, remove the pair where we've excluded the larger barcode
    for minor_barcode in barcode_multiplets.keys():
        if minor_barcode not in barcode_multiplets:
            # Because we've popped it off before we got here
            continue
        major_barcode = barcode_multiplets[minor_barcode]
        if major_barcode in barcode_multiplets and barcode_multiplets[major_barcode] == minor_barcode:
            if barcode_counts[major_barcode] > barcode_counts[minor_barcode]:
                barcode_multiplets.pop(major_barcode)
            else:
                barcode_multiplets.pop(minor_barcode)

    # Post-screen barcode multiplets for those where the major barcode is itself
    # linked to another barcode
    for minor_barcode, major_barcode in barcode_multiplets.iteritems():
        if major_barcode in barcode_multiplets:
            major_barcode = barcode_multiplets[major_barcode]
            barcode_multiplets[minor_barcode] = major_barcode

    # Generate the exclusions.  Note we write it out once per species since
    # cell calling is species-specific but these exclusions are not.
    ref = ReferenceManager(args.reference_path)
    species_list = ref.list_species()
    excluded_barcodes = {
        "label": "whitelist_contam",
        "data": {species: barcode_multiplets for species in species_list}
    }
    with open(outs.barcode_multiplets, "w") as outfile:
        outfile.write(json.dumps(excluded_barcodes))

    # Generate some reporting metrics
    summary_metrics = {
        "putative_barcode_multiplets_found": len(barcode_multiplets),
    }
    with open(outs.barcode_multiplets_summary, "w") as outfile:
        outfile.write(json.dumps(summary_metrics))


def main(args, outs):
    with open(args.barcode_counts, "r") as infile:
        barcode_counts = Counter(json.load(infile))

    valid_barcodes = barcode_counts.keys()
    part_a_seqs, part_c_seqs, part_b_seqs, gem_group_seqs = query_barcode_subsequences(
        valid_barcodes
    )

    part_a_count = max([len(part_a_seqs[c]) for c in part_c_seqs])
    part_b_count = max([len(part_b_seqs[c]) for c in part_c_seqs])
    part_c_count = len(part_c_seqs)

    index_by_part_a = {
        part_c: {
            part_a: i for i, part_a in enumerate(part_a_seqs[part_c])
        } for part_c in part_c_seqs
    }
    index_by_part_b = {
        part_c: {
            part_b: i for i, part_b in enumerate(part_b_seqs[part_c])
        } for part_c in part_c_seqs
    }
    index_by_part_c = {
        part_c: i for i, part_c in enumerate(part_c_seqs)
    }

    part_a_linkage_matrix = np.zeros((part_c_count, part_b_count, part_a_count, part_a_count), dtype=np.uint32)
    part_b_linkage_matrix = np.zeros((part_c_count, part_a_count, part_b_count, part_b_count), dtype=np.uint32)

    stops = {}

    for fragment_list in grouped_fragments_from_contig(
        args.contig, args.fragments, args.fragments_index
    ):
        for _, start, stop, barcode, _ in fragment_list:
            part_a, part_c, part_b, gem_group = split_barcode(barcode, return_gg=True)
            if gem_group != args.gem_group:
                # Only include fragments with the required input gem group.
                continue

            part_a_index = index_by_part_a[part_c][part_a]
            part_b_index = index_by_part_b[part_c][part_b]
            part_c_index = index_by_part_c[part_c]

            # Periodically clean up old keys to reduce memory footprint
            if len(stops) > MAXIMUM_POSITION_SIZE:
                for key in stops.keys():
                    if key < (start - MAXIMUM_FRAGMENT_SIZE):
                        stops.pop(key)

            if stop not in stops:
                stops[stop] = Counter()
            stops[stop][(part_c_index, part_a_index, part_b_index)] += 1

            if start in stops and (len(stops[stop]) + len(fragment_list)) <= MAXIMUM_PILEUP:
                for (last_part_c_index, last_part_a_index, last_part_b_index), count in stops[start].iteritems():
                    if part_c_index != last_part_c_index:
                        continue

                    if part_b_index == last_part_b_index:
                        part_a_linkage_matrix[part_c_index, part_b_index, part_a_index, last_part_a_index] += 1
                        if part_a_index != last_part_a_index:
                            part_a_linkage_matrix[part_c_index, part_b_index, last_part_a_index, part_a_index] += 1

                    if part_a_index == last_part_a_index:
                        part_b_linkage_matrix[part_c_index, part_a_index, part_b_index, last_part_b_index] += 1
                        if part_b_index != last_part_b_index:
                            part_b_linkage_matrix[part_c_index, part_a_index, last_part_b_index, part_b_index] += 1

    outfile = gzip.GzipFile(outs.part_a_linkage_matrix, "w")
    np.save(outfile, part_a_linkage_matrix)
    outfile.close()

    outfile = gzip.GzipFile(outs.part_b_linkage_matrix, "w")
    np.save(outfile, part_b_linkage_matrix)
    outfile.close()
