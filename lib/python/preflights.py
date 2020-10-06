import socket
import os
import numpy as np
import pandas as pd
import re

import martian
import tenkit.reference as tk_ref
import tenkit.preflight as tk_preflight
from constants import ALLOWED_FACTORIZATIONS, FRAGMENTS_SCAN_SIZE, NO_BARCODE
from utils import motif_format_checker
from tools.regions import bed_format_checker, is_overlapping
from tools import ReferenceManager
from tools.io import open_fragment_file, parse_aggr_csv

def check_force_cells(force_cells, ulimit=20000):
    """check if force cells is correctly specified"""
    if force_cells is not None:
        if len(force_cells) == 0:
            martian.exit("force_cells must be a non-empty dictionary.")
        for force_cells_k in force_cells.keys():
            if (force_cells[force_cells_k] < 1 or force_cells[force_cells_k] > ulimit):
                martian.exit("MRO parameter force-cells[{}] must be a positive integer <= 20000: {}".
                             format(force_cells_k, force_cells[force_cells_k]))


def check_sample_id(sample_id):
    """check sample id is correctly formatted"""
    if sample_id is not None:
        if not re.match("^[\w-]+$", sample_id):
            martian.exit("Sample name may only contain letters, numbers, underscores, and dashes: {}".format(sample_id))

def exists_and_readable(fname, key, isfile=True):
    """check if file exists and is readable"""
    if not os.path.exists(fname):
        martian.exit("Specified {} file: {} does not exist".format(key, fname))
    if isfile:
        if not os.path.isfile(fname):
            martian.exit("Specified {} file: {} is not a valid file".format(key, fname))
    if not os.access(fname, os.R_OK):
        martian.exit("{} is not readable, please check file permissions".format(fname))


def check_singlecell_format(fname, species_list, whitelist=None, allow_multi_gem_groups=True):
    """check if the singlecell.csv file has at minimum:
    the barcodes column with barcodes on the whitelist
    and is_{}_cell_barcode columns with right value type and ranges"""
    exists_and_readable(fname, "singlecell")
    scdf = pd.read_csv(fname, sep=',')
    if scdf.isnull().values.any():
        martian.exit("Nans are not acceptable values in {}. This file is likely ill-formatted, with missing or empty entries in some columns for some row(s)".format(fname))

    # expect barcodes
    if "barcode" not in scdf.columns:
        martian.exit("barcode is a required column in {}".format(fname))
    else:
        if scdf['barcode'].dtype != 'O':
            martian.exit("barcode column has keys of unexpected type in {}".format(fname))

    # barcodes should be on whitelist
    if whitelist is not None or not allow_multi_gem_groups:
        if whitelist is not None:
            assert type(whitelist) == set
        gemgroups = set()

        for barcode in scdf['barcode']:
            if barcode != NO_BARCODE:
                bc, gemgroup = barcode.split("-")
            else:
                bc = barcode
            if whitelist is not None:
                if barcode != NO_BARCODE and bc not in whitelist:
                    martian.exit("Found barcode {} in {} that is not on the barcode whitelist".format(barcode, fname))
            if not allow_multi_gem_groups:
                if barcode != NO_BARCODE:
                    gemgroups.add(gemgroup)
                    if len(gemgroups) > 1:
                        martian.exit("Found barcodes from multiple gem groups in {}. This data is not currently supported on the pipeline".format(fname))

    # expect cell annotation columns to match reference
    for species in species_list:
        species_cells = "is_{}_cell_barcode".format(species)
        if species_cells not in scdf.columns:
            msg = "expecting cell barcodes of species {} in {}.".format(species, fname) if len(species_list) > 1 else "expecting 'is__cell_barcode' column in {}".format(fname)
            martian.exit(msg + " This could be because of specifying the wrong reference.")
        else:
            if scdf[species_cells].dtype != 'int64':
                martian.exit("{} column has keys of unexpected type in {}".format(species_cells, fname))
            if np.sum(scdf[species_cells] > 1) > 0:
                martian.exit("{} column has keys of unexpected values in {}. Allowed values are 0 or 1".format(species_cells, fname))

def check_aggr_csv(aggr_csv, reference_path, cursory=False):
    """Check aggr csv has correct columns, then progressively stronger checks on duplicates and formating of files.
    These stronger checks are enabled by default, unless you want to test the basic minimum, for example in reanalyzer"""
    contig_manager = ReferenceManager(reference_path)

    # aggr_csv checks
    exists_and_readable(aggr_csv, "aggr_csv")

    if cursory:
        nlibs, library_info, msg = parse_aggr_csv(aggr_csv, whitelist=["library_id"], blacklist=None)
    else:
        nlibs, library_info, msg = parse_aggr_csv(aggr_csv)
    if msg is not None:
        martian.exit(msg)

    # At least one library should be there
    if nlibs == 0:
        martian.exit("aggregation csv does not include any library. Provide at least two libraries.")

    if cursory:
        return
    # Enable aggr(count1) to run
    if nlibs == 1:
        martian.log_info("Aggregator should be run on more than one library")

    # avoid aggr of duplicate files (assessed by filename).
    species_list = contig_manager.list_species()
    for aggr_key in library_info[1]:  # at least one library is present
        files = {}
        for lib_id in library_info:
            fname = library_info[lib_id][aggr_key]
            if fname in files:
                martian.exit("File {} already specified for a different library under {}".format(fname, aggr_key))

            # singlecell.csv should contain 'barcode' and 'is_{}_cell_barcode' columns with the correct type
            if aggr_key == "cells":
                check_singlecell_format(fname, species_list, allow_multi_gem_groups=False)

            # peaks.bed need to be formatted correctly with right contigs if provided in aggr.csv
            # also check if peaks are non overlapping
            if aggr_key == "peaks":
                exists_and_readable(fname, "peaks")
                bed_format_checker(fname, contig_manager.fasta_index)
                contain_three_columns(fname)
                if is_overlapping(fname):
                    martian.exit("{} contains overlapping peak regions".format(fname))

            # checks on fragments
            contig_lens = contig_manager.get_contig_lengths()
            if aggr_key == "fragments":
                observed_gem_groups = set()
                observed_species = set()
                exists_and_readable(fname, "fragments")
                en = 0
                for chrom, start, stop, bc, _ in open_fragment_file(fname):
                    if en >= FRAGMENTS_SCAN_SIZE:
                        break
                    spec = chrom.split("_")
                    observed_species.add(spec[0] if spec[0] != chrom else "")
                    observed_gem_groups.add(bc.split("-")[1])
                    if chrom not in contig_lens:
                        martian.exit("fragment {}:{}-{} in {} is mapped to a contig not in the reference".format(chrom, start, stop, fname))
                    if stop > contig_lens[chrom]:
                        martian.exit("fragment {}:{}-{} boundaries exceed contig size ({} bp)".format(chrom, start, stop, contig_lens[chrom]))
                    en += 1
                for species in observed_species:
                    if species not in species_list:
                        martian.exit("{} contains fragments mapped to species not recognized in the reference".format(fname))
                if len(observed_gem_groups) > 1:
                    martian.exit("multiple gem groups present in {}, likely generated in a previous aggregation run".format(fname))


def check_filehandle_limit():
    """checks file handles"""
    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

def check_factorization(factorization):
    """checks if factorization is valid"""
    if factorization is not None:
        if len(factorization) == 0:
            martian.exit("factorization must be a non-empty list.")
        if not all(elem in ALLOWED_FACTORIZATIONS for elem in factorization):
            martian.exit("Unsupported factorization provided. Options are {}.".
                         format(", ".join(ALLOWED_FACTORIZATIONS)))


def check_vmem_for_reference(ref):
    """Finds out vmem required to load a genome reference fully"""
    refpath = os.path.join(ref, "fasta", "genome.fa")
    if not os.path.exists(refpath):
        hostname = socket.gethostname()
        martian.exit("Your reference does not contain the expected files, or they are not readable. Please check your reference folder on {}.".format(hostname))
    refsize = os.path.getsize(refpath) / 1e9
    vmem_gb = int(np.ceil(refsize)) + 4
    return vmem_gb

def check_reference_format(reference_path):
    """Check file formats for files present in the reference"""
    try:
        contig_manager = ReferenceManager(reference_path)
    except Exception as e:
        martian.exit("Contig manager could not be initialized, Error:\n%s" % str(e))

    # formatting
    error_msg = contig_manager.verify_contig_defs()
    if error_msg is not None:
        martian.exit(error_msg)

    # filecheck
    contig_manager.genes

    # check if motif file is in right format (naming convention)
    if len(contig_manager.list_species()) == 1:
        motif_format_checker(contig_manager.motifs)

    # checks for valid bed file formats in regions/
    faidx_file = os.path.join(reference_path, 'fasta', 'genome.fa.fai')

    bed_format_checker(contig_manager.tss_track, faidx_file)
    bed_format_checker(contig_manager.transcripts_track, faidx_file)
    bed_format_checker(contig_manager.ctcf_track, faidx_file)
    bed_format_checker(contig_manager.blacklist_track, faidx_file)
    bed_format_checker(contig_manager.dnase_track, faidx_file)
    bed_format_checker(contig_manager.enhancer_track, faidx_file)
    bed_format_checker(contig_manager.promoter_track, faidx_file)


def check_refdata(reference_path, max_contigs=None):
    """Does various checks on reference files structure"""
    hostname = socket.gethostname()

    def pathjoin(relpath):
        return os.path.join(reference_path, relpath)

    def pathexists(relpath):
        return os.path.exists(pathjoin(relpath))

    def getmodtime(relpath):
        return os.path.getmtime(pathjoin(relpath))

    # Determine if the reference package is a known 10X reference package
    genome = tk_ref.get_genome(reference_path)

    if genome is not None and tk_ref.is_tenx(reference_path):
        if not pathexists("version"):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

    # Known genomes get a more stringent check
    if not pathexists("fasta/") or not pathexists("genes/"):
        return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

    if not pathexists("fasta/genome.fa.flat") or not pathexists("fasta/genome.fa.gdx"):
        return False, "Your reference doesn't appear to be indexed. Please run the mkref tool."

    if getmodtime("fasta/genome.fa.flat") < getmodtime("fasta/genome.fa") or getmodtime("fasta/genome.fa.gdx") < getmodtime("fasta/genome.fa"):
        return False, "Timestamps on the FASTA index files have changed. This may have happened during a copy or move operation on the reference files. Please reindex your reference if needed, or modify their timestamp using touch linux command"

    fasta = tk_ref.open_reference(reference_path)
    num_contigs = len(fasta.keys())

    if max_contigs is not None and num_contigs > max_contigs:
        return False, "Long Ranger supports a maximum of %d reference contigs. Your reference contains %d. Please combine small contigs into a larger contig separated by N's." % (max_contigs, num_contigs)

    max_len = max([len(v) for (k, v) in fasta.iteritems()])

    logging = "reference path %s on %s contains genome: %s." % (reference_path, hostname, str(genome))
    logging += "reference contains %d contigs. max contig length: %d." % (num_contigs, max_len)

    if max_len >= (1 << 29):
        return False, "Reference contains a contig longer than 536.8Mb (2^29 bp), which is not supported due to limitations of the .bai format. Please split this contig."

    # Check for ":" in contig names -- this will break things horribly
    has_colons = any(":" in ctg_name for ctg_name in fasta.keys())
    if has_colons:
        return False, "Reference names contain colon characters: ':'. References names containing colons are not supported."

    return True, logging

def contain_three_columns(in_file, check_top_n = 100):

    with open(in_file) as bediter:
        for i, row in enumerate(bediter):
            fields = row.strip().split('\t')
            if len(fields) != 3:
                martian.exit("Peak input file must contain only 3 columns")

            if i > check_top_n:
                break
        return
