#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

import os

import martian
import tenkit.preflight as tk_preflight
from tools import ReferenceManager
from preflights import (check_refdata, check_vmem_for_reference, check_filehandle_limit,
                        check_sample_id, check_reference_format, check_aggr_csv,
                        check_force_cells, check_singlecell_format, exists_and_readable, contain_three_columns)
from tools.regions import bed_format_checker, is_overlapping
from tools.io import open_fragment_file, parsed_fragments_from_contig
from barcodes import load_barcode_whitelist
from constants import FRAGMENTS_SCAN_SIZE

__MRO__ = """
stage ATAC_REANALYZER_PREFLIGHT(
    in  string     sample_id,
    in  string     reference_path,
    in  string     barcode_whitelist,
    in  bed        peaks,
    in  csv        parameters,
    in  map        force_cells,
    in  csv        cell_barcodes,
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    in  csv        aggregation_csv,
    in  bool       check_executables,
    src py         "stages/preflight/atac_reanalyzer",
) split (
)
"""

def split(args):
    vmem_gb = check_vmem_for_reference(args.reference_path)
    return {'chunks': [], 'join': {'__vmem_gb': vmem_gb}}

def main(args, outs):
    martian.throw("main not implemented")

def join(args, outs, chunk_defs, chunk_outs):
    # Sample ID / pipestance name
    check_sample_id(args.sample_id)

    # force_cells
    check_force_cells(args.force_cells, ulimit=10000000)  # allow arbitrarily large limit for reanalyzer

    # # Reference
    # ref directory structure and timestamps
    ok, msg = check_refdata(args.reference_path, max_contigs=None)
    if ok:
        martian.log_info(msg)
    else:
        martian.exit(msg)

    # formatting
    check_reference_format(args.reference_path)
    contig_manager = ReferenceManager(args.reference_path)

    # peaks format check and nonoverlapping
    if args.peaks is None:
        martian.exit("peaks file not provided")
    exists_and_readable(args.peaks, "peaks")
    bed_format_checker(args.peaks, contig_manager.fasta_index)
    contain_three_columns(args.peaks)
    if is_overlapping(args.peaks):
        martian.exit("{} contains overlapping peak regions".format(args.peaks))

    # check parameters files
    if args.parameters is not None:
        if not os.path.exists(args.parameters):
            martian.exit("{} does not exist".format(args.parameters))

    # fragments checks
    whitelist_barcodes = load_barcode_whitelist(args.barcode_whitelist)
    species_list = contig_manager.list_species()
    observed_gem_groups = set()
    observed_species = set()
    if args.fragments is None:
        martian.exit("fragments file not provided")
    exists_and_readable(args.fragments, "fragments")
    contig_lens = contig_manager.get_contig_lengths()
    # check bounds and matching contigs in reference and species
    for chrom, start, stop, bc, _ in open_fragment_file(args.fragments):
        spec = chrom.split("_")
        observed_species.add(spec[0] if spec[0] != chrom else "")
        barcode, gem_group = bc.split("-")
        observed_gem_groups.add(gem_group)
        if args.check_executables:  # run this only non-locally
            if barcode not in whitelist_barcodes:
                martian.exit("{} is not a valid whitelist barcode".format(barcode))
            if chrom not in contig_lens:
                martian.exit("contig {} not present in reference".format(chrom))
            if stop > contig_lens[chrom]:
                martian.exit("fragment {}:{}-{} boundaries exceed contig size ({} bp)".format(chrom, start, stop, contig_lens[chrom]))
    # ensure fragments are on the correct reference
    for species in observed_species:
        if species not in species_list:
            martian.exit("{} contains fragments mapped to species not recognized in the reference".format(args.fragments))
    if len(observed_gem_groups) > 1:
        martian.log_info("multiple gem groups present in {}, likely generated in a previous aggregation run".format(args.fragments))

    # fragments index is synced with fragments
    if args.fragments_index is None:
        martian.exit("fragments index file not provided")
    if not os.path.exists(args.fragments_index):
        martian.exit("{} does not exist".format(args.fragments_index))
    try:
        all_contigs = contig_manager.primary_contigs(allow_sex_chromosomes=True)
        for contig in all_contigs:
            en = 0
            for chrom, start, end, bc, dups in parsed_fragments_from_contig(contig, args.fragments, index=args.fragments_index):
                if en >= FRAGMENTS_SCAN_SIZE:
                    break
                en += 1
    except:
        martian.exit("fragments index is not in sync with the fragments file")

    # aggr csv checks
    if args.aggregation_csv is not None:
        check_aggr_csv(args.aggregation_csv, args.reference_path, cursory=True)

    # cell barcode checks
    if args.cell_barcodes is not None:
        if not os.path.exists(args.cell_barcodes):
            martian.exit("{} does not exist".format(args.cell_barcodes))
        check_singlecell_format(args.cell_barcodes, species_list, whitelist_barcodes)

    # Open file handles limit
    if args.check_executables:
        check_filehandle_limit()

    martian.log_info(tk_preflight.record_package_versions())
