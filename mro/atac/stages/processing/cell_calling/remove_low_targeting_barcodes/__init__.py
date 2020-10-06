"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Identifies barcodes with low targeting (highly likely to be from dead cells), with no
limitation to accessible regions.

Outputs a json of species-specific barcodes to be excluded from later cell calling.
"""
from __future__ import division

import json
import pickle
from collections import Counter

from tools import ReferenceManager, parsed_fragments_from_contig
from tools.regions import get_target_regions, fragment_overlaps_target

__MRO__ = """
stage REMOVE_LOW_TARGETING_BARCODES(
    in  bed        peaks,
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  string     reference_path,
    out json       barcode_counts,
    out json       low_targeting_barcodes,
    out json       low_targeting_summary,
    out json       fragment_lengths,
    out json       covered_bases,
    src py         "stages/processing/cell_calling/remove_low_targeting_barcodes",
) split (
    in  string     contig,
    out pickle     fragment_counts,
    out pickle     targeted_counts,
    out int        peak_coverage,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""

DISTANCE = 2000
PADDING_VALUES = [0, 50000]


def count_covered_bases(start_iter, stop_iter, contig_len, padding=1000):
    """
    Loop through an input file, counting up the per-barcode total coverage and maximum
    coverage
    """
    occupancy = 0
    last_start = None
    last_stop = None

    for start, stop in zip(start_iter, stop_iter):
        start = max(0, start - padding)
        stop = min(contig_len, stop + padding)

        if last_start is None:
            last_start = start
            last_stop = stop
        else:
            if start < last_stop:
                last_stop = max(stop, last_stop)
            else:
                occupancy += last_stop - last_start
                last_start = start
                last_stop = stop

    if last_start is not None:
        occupancy += last_stop - last_start

    return occupancy


def split(args):
    if args.fragments is None:
        return {"chunks": [], "join": {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    chunks = []
    for contig in all_contigs:
        chunks.append({"contig": contig, "__mem_gb": 5})

    return {"chunks": chunks, "join": {"__mem_gb": 5}}


def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()
    if args.fragments is None:
        outs.low_targeting_barcodes = None
        outs.low_targeting_summary = None
        return

    # Merge the chunk inputs
    ref = ReferenceManager(args.reference_path)
    species_list = ref.list_species()

    barcode_counts_by_species = {species: Counter() for species in species_list}
    targeted_counts_by_species = {species: Counter() for species in species_list}

    peak_bp_by_species = {species: 0 for species in species_list}
    genome_bp_by_species = {species: 0 for species in species_list}

    fragment_lengths = {padding: Counter() for padding in PADDING_VALUES}
    covered_bases = {padding: Counter() for padding in PADDING_VALUES}

    for chunk_in, chunk_out in zip(chunk_defs, chunk_outs):
        species = ref.species_from_contig(chunk_in.contig)

        with open(chunk_out.fragment_counts, "r") as infile:
            barcode_counts_by_species[species] += pickle.load(infile)
        with open(chunk_out.targeted_counts, "r") as infile:
            targeted_counts_by_species[species] += pickle.load(infile)

        with open(chunk_out.fragment_lengths, "r") as infile:
            data = pickle.load(infile)
            for padding in PADDING_VALUES:
                fragment_lengths[padding] += data[padding]

        with open(chunk_out.covered_bases, "r") as infile:
            data = pickle.load(infile)
            for padding in PADDING_VALUES:
                covered_bases[padding] += data[padding]

        peak_bp_by_species[species] += chunk_out.peak_coverage
        genome_bp_by_species[species] += ref.contig_lengths[chunk_in.contig]

    frac_genome_in_peaks_by_species = {
        species: peak_bp_by_species[species] / genome_bp_by_species[species]
        for species in species_list
    }

    # Identify barcodes that have lower fraction of reads overlapping peaks than the
    # genomic coverage of the peaks
    low_targeting_barcodes = {
        "label": "low_targeting",
        "data": {species: {} for species in species_list}
    }
    for species in species_list:
        for barcode, total_count in barcode_counts_by_species[species].iteritems():
            barcode_frac_peaks = (
                targeted_counts_by_species[species][barcode] / total_count
            )
            if barcode_frac_peaks < frac_genome_in_peaks_by_species[species]:
                low_targeting_barcodes["data"][species][barcode] = barcode_frac_peaks

    # Sum up the total fragment counts per barcode across all species
    total_barcode_counts = Counter()
    for species, barcode_counts in barcode_counts_by_species.iteritems():
        total_barcode_counts += barcode_counts
    with open(outs.barcode_counts, "w") as outfile:
        outfile.write(json.dumps(total_barcode_counts, indent=4))

    summary_data = {}
    for species in species_list:
        key_suffix = "" if len(species_list) == 1 else "_{}".format(species)
        summary_data["number_of_low_targeting_barcodes{}".format(key_suffix)] = len(
            low_targeting_barcodes["data"][species]
        )
        summary_data[
            "fraction_of_genome_within_{}bp_of_peaks{}".format(DISTANCE, key_suffix)
        ] = frac_genome_in_peaks_by_species[species]
    with open(outs.low_targeting_summary, "w") as outfile:
        outfile.write(json.dumps(summary_data, indent=4))
    with open(outs.low_targeting_barcodes, "w") as outfile:
        outfile.write(json.dumps(low_targeting_barcodes, indent=4))
    with open(outs.fragment_lengths, "w") as outfile:
        outfile.write(json.dumps(fragment_lengths, indent=4))
    with open(outs.covered_bases, "w") as outfile:
        outfile.write(json.dumps(covered_bases, indent=4))


def main(args, outs):
    ref = ReferenceManager(args.reference_path)
    contig_len = ref.get_contig_lengths()[args.contig]

    with open(args.peaks, "r") as infile:
        peak_regions = get_target_regions(infile)

    fragment_counts = Counter()
    targeted_counts = Counter()

    cumulative_fragment_length = {padding: Counter() for padding in PADDING_VALUES}
    covered_bases = {padding: {} for padding in PADDING_VALUES}

    for contig, start, stop, barcode, _ in parsed_fragments_from_contig(
        args.contig, args.fragments, args.fragments_index
    ):
        fragment_counts[barcode] += 1
        if fragment_overlaps_target(contig, start, stop, peak_regions):
            targeted_counts[barcode] += 1

        for padding in PADDING_VALUES:
            adj_start = max(0, start - padding)
            adj_stop = min(contig_len, stop + padding)

            cumulative_fragment_length[padding][barcode] += adj_stop - adj_start

            if barcode not in covered_bases[padding]:
                # Total # of covered bases, current start/stop of active region
                covered_bases[padding][barcode] = 0, None, None

            current_covered, active_start, active_stop = covered_bases[padding][barcode]
            if active_start is None:
                active_start = adj_start
                active_stop = adj_stop
            else:
                if adj_start < active_stop:
                    active_stop = max(adj_stop, active_stop)
                else:
                    current_covered += active_stop - active_start
                    active_start = adj_start
                    active_stop = adj_stop
            covered_bases[padding][barcode] = current_covered, active_start, active_stop

    final_covered = {padding: Counter() for padding in PADDING_VALUES}
    for padding in PADDING_VALUES:
        for barcode in covered_bases[padding]:
            current_covered, active_start, active_stop = covered_bases[padding][barcode]
            if active_start is None:
                final_covered[padding][barcode] = current_covered
            else:
                final_covered[padding][barcode] = current_covered + active_stop - active_start

    with open(outs.fragment_counts, "w") as outfile:
        pickle.dump(fragment_counts, outfile)
    with open(outs.targeted_counts, "w") as outfile:
        pickle.dump(targeted_counts, outfile)
    with open(outs.fragment_lengths, "w") as outfile:
        pickle.dump(cumulative_fragment_length, outfile)
    with open(outs.covered_bases, "w") as outfile:
        pickle.dump(final_covered, outfile)

    outs.peak_coverage = 0
    if args.contig in peak_regions:
        contig_len = ref.contig_lengths[args.contig]
        outs.peak_coverage = count_covered_bases(
            peak_regions[args.contig].starts,
            peak_regions[args.contig].ends,
            contig_len,
            padding=DISTANCE,
        )
