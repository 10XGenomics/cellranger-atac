"""
Generate a connectivity matrix between barcodes with some minimum number of fragments,
connecting barcodes when two fragments from the barcodes are directly adjacent,
indicative of a shared cut site.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
from __future__ import absolute_import, division, print_function

from tools import ReferenceManager, grouped_fragments_from_contig
import numpy as np
import json
from collections import Counter
from barcodes import (get_barcode_gem_group, query_barcodes_and_gem_groups,
                      split_barcode_and_gem_group, merge_barcode_and_gem_group)
import martian

__MRO__ = """
stage REMOVE_GEL_BEAD_DOUBLET_BARCODES(
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  string     reference_path,
    in  json       barcode_counts,
    out json       gel_bead_doublet_barcodes,
    out json       gel_bead_doublet_summary,
    out csv        connect_matrix,
    src py         "stages/processing/cell_calling/remove_gel_bead_doublet_barcodes",
) split (
    in  string     contig,
    in  file       valid_barcodes,
) using (
    mem_gb = 4,
    volatile = strict,
)
"""

# Minimum number of fragments for a barcode to be included in the connectivity matrix
MINIMUM_COUNTS = 250
MAXIMUM_BARCODES = 15000

MAXIMUM_FRAGMENT_SIZE = 800
MAXIMUM_POSITION_SIZE = 2500

MAXIMUM_PILEUP = 20


def split(args):
    if args.fragments is None:
        return {"chunks": [], "join": {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    with open(args.barcode_counts, "r") as infile:
        barcode_counts = Counter(json.load(infile))
    barcode_array = np.array([bc for bc in barcode_counts])
    gem_group_array = np.array([get_barcode_gem_group(bc) for bc in barcode_counts])
    gem_groups = set(gem_group_array)
    frag_count_array = np.array([barcode_counts[bc] for bc in barcode_array])

    valid_barcodes = list()
    for gem_group in gem_groups:
        count_mask = (frag_count_array > MINIMUM_COUNTS) & (gem_group_array == gem_group)
        # find at most top N barcodes
        topN_indices = barcode_array[count_mask].argsort()[-min(MAXIMUM_BARCODES, len(count_mask)):]
        valid_barcodes.extend(list(barcode_array[count_mask][topN_indices]))

    # mem allocs
    JOIN_LOAD_FACTOR = 2
    BUFFER_GB = 2
    BYTES_PER_ENTRY = 4  # this depends on the dtype
    chunk_mem_gb = BUFFER_GB + np.ceil(BYTES_PER_ENTRY * len(gem_groups) * MAXIMUM_BARCODES ** 2 / 1024 ** 3).astype('int32')
    join_mem_gb = BUFFER_GB + np.ceil(JOIN_LOAD_FACTOR * BYTES_PER_ENTRY * len(gem_groups) * MAXIMUM_BARCODES ** 2 / 1024 ** 3).astype('int32')

    valid_barcodes_path = martian.make_path("valid_barcodes.txt")
    with open(valid_barcodes_path, 'w') as f:
        f.write(",".join(valid_barcodes))

    chunks = []
    for contig in all_contigs:
        chunks.append(
            {
                "contig": contig,
                "valid_barcodes": valid_barcodes_path,
                "__mem_gb": chunk_mem_gb,
            }
        )

    return {"chunks": chunks, "join": {"__mem_gb": join_mem_gb}}


def join(args, outs, chunk_defs, chunk_outs):
    if args.fragments is None:
        outs.connect_matrix = None
        outs.gel_bead_doublet_summary = None
        outs.gel_bead_doublet_barcodes = None
        return

    with open(args.barcode_counts, "r") as infile:
        barcode_counts = Counter(json.load(infile))

    with open(chunk_defs[0].valid_barcodes, 'r') as f:
        valid_barcodes = np.array(f.readlines()[0].strip("\n").split(","))

    barcode_seqs, gem_groups = query_barcodes_and_gem_groups(valid_barcodes)
    barcode_seq_count = max([len(barcode_seqs[gg]) for gg in gem_groups])
    n_gem_groups = len(gem_groups)
    index_by_barcode = {gg: {bc: i for i, bc in enumerate(barcode_seqs[gg])} for gg in gem_groups}
    index_by_gg = {gg: i for i, gg in enumerate(gem_groups)}

    connect_matrix = np.zeros(
        (n_gem_groups, barcode_seq_count, barcode_seq_count), dtype=np.uint32
    )

    # This can be memory intensive due to loading the same amount of memory
    for chunk_out in chunk_outs:
        with open(chunk_out.connect_matrix, "r") as infile:
            connect_matrix += np.load(infile)

    # Write out the raw matrix
    with open(outs.connect_matrix, "w") as outfile:
        for gg in gem_groups:
            outfile.write(",".join([merge_barcode_and_gem_group(bc, gg) for bc in barcode_seqs[gg]]))
            outfile.write("\n")
            for i in range(len(barcode_seqs[gg])):
                outfile.write(",".join((str(count) for count in connect_matrix[index_by_gg[gg], i, :])))
                outfile.write("\n")

    # Identify mutual nearest neighbors as putative doublets
    putative_doublets = []
    for barcode in valid_barcodes:
        bc_seq, gg = split_barcode_and_gem_group(barcode)
        gg_index = index_by_gg[gg]
        bc_index = index_by_barcode[gg][bc_seq]
        neighbor = nearest_neighbor(connect_matrix, bc_index, gg_index)
        if nearest_neighbor(connect_matrix, neighbor, gg_index) == bc_index:
            if bc_index < neighbor:
                putative_doublets.append((barcode,
                                          merge_barcode_and_gem_group(barcode_seqs[gg][neighbor], gg)))

    # Generate the exclusions.  Note we write it out once per species since
    # cell calling is species-specific but these exclusions are not.
    ref = ReferenceManager(args.reference_path)
    species_list = ref.list_species()
    excluded_barcodes = {
        "label": "gel_bead_doublet",
        "data": {species: {} for species in species_list}
    }
    for pair in putative_doublets:
        if barcode_counts[pair[0]] < barcode_counts[pair[1]]:
            excluded_bc, major_bc = pair
        else:
            major_bc, excluded_bc = pair
        for species in species_list:
            excluded_barcodes["data"][species][excluded_bc] = major_bc
    with open(outs.gel_bead_doublet_barcodes, "w") as outfile:
        outfile.write(json.dumps(excluded_barcodes))

    estimated_doublet_gelbeads = len(putative_doublets)

    metrics = {"putative_gelbead_doublets_found": estimated_doublet_gelbeads}
    with open(outs.gel_bead_doublet_summary, "w") as outfile:
        outfile.write(json.dumps(metrics))


def nearest_neighbor(pair_matrix, bc_index, gg_index):
    """Find the non-self index in a connectivity matrix with the highest value.
    Do this in a gem group specific way"""
    row = pair_matrix[gg_index, bc_index, :]
    indices = np.arange(pair_matrix.shape[1])
    mask = indices != bc_index
    return indices[mask][np.argmax(row[mask])]


def main(args, outs):
    with open(args.valid_barcodes, 'r') as f:
        valid_barcodes = f.readlines()[0].strip("\n").split(",")

    barcode_seqs, gem_groups = query_barcodes_and_gem_groups(valid_barcodes)
    barcode_seq_count = max([len(barcode_seqs[gg]) for gg in gem_groups])
    n_gem_groups = len(gem_groups)
    index_by_barcode = {gg: {bc: i for i, bc in enumerate(barcode_seqs[gg])} for gg in gem_groups}
    index_by_gg = {gg: i for i, gg in enumerate(gem_groups)}

    connect_matrix = np.zeros(
        (n_gem_groups, barcode_seq_count, barcode_seq_count), dtype=np.uint32
    )
    stops = {}

    for fragment_list in grouped_fragments_from_contig(
        args.contig, args.fragments, args.fragments_index
    ):
        for _, start, stop, barcode, _ in fragment_list:
            barcode_seq, gem_group = split_barcode_and_gem_group(barcode)
            if gem_group not in index_by_gg:
                continue
            if barcode_seq not in index_by_barcode[gem_group]:
                continue

            # Periodically clean up old keys to reduce memory footprint
            if len(stops) > MAXIMUM_POSITION_SIZE:
                for key in stops.keys():
                    if key < (start - MAXIMUM_FRAGMENT_SIZE):
                        stops.pop(key)

            index1 = index_by_barcode[gem_group][barcode_seq]
            gg_index = index_by_gg[gem_group]
            if stop not in stops:
                stops[stop] = Counter()
            stops[stop][barcode] += 1

            if start in stops and (len(stops[stop]) + len(fragment_list)) <= MAXIMUM_PILEUP:
                for bc, count in stops[start].iteritems():
                    bc_seq, other_gem_group = split_barcode_and_gem_group(bc)
                    if gem_group != other_gem_group:
                        continue

                    index2 = index_by_barcode[gem_group][bc_seq]
                    connect_matrix[gg_index, index1, index2] += count
                    if index1 != index2:
                        connect_matrix[gg_index, index2, index1] += count

    with open(outs.connect_matrix, "w") as outfile:
        np.save(outfile, connect_matrix)
