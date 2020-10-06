#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

import os
import re
import socket

import martian
import tenkit.fasta as tk_fasta
import tenkit.preflight as tk_preflight
from preflights import (check_refdata, check_vmem_for_reference, check_filehandle_limit,
                        check_sample_id, check_force_cells, check_factorization,
                        check_reference_format)
from constants import BCL2FASTQ_SEQNAMES

__MRO__ = """
stage ATAC_COUNTER_PREFLIGHT(
    in  string   sample_id,
    in  string   fastq_mode,
    in  map[]    sample_def,
    in  string   reference_path,
    in  map      force_cells,
    in  string[] factorization,
    in  map      downsample,
    in  bool     check_executables,
    in  map      trim_def,
    src py       "stages/preflight/atac_counter",
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
    check_force_cells(args.force_cells)

    # downsample
    if args.downsample is not None:
        if len(args.downsample) == 0:
            martian.exit("downsample must be a non-empty dictionary.")
        keys = args.downsample.keys()
        if len(keys) > 1:
            martian.exit("Please supply either subsample_rate or gigabases but not both.")
        key = keys[0]
        if not (key in ['subsample_rate', 'gigabases']):
            martian.exit("Please supply either subsample_rate or gigabases as the downsample argument. '%s' is invalid" % key)
        value = args.downsample[key]
        bad_value = False
        try:
            float(value)
            bad_value = value < 1e-12
        except ValueError:
            bad_value = True
        if bad_value:
            martian.exit("Command line argument for downsampling must be a positive number")

    # FASTQ mode
    if args.fastq_mode is not None:
        if args.fastq_mode not in ['ILMN_BCL2FASTQ', 'BCL_PROCESSOR']:
            martian.exit("Unsupported fastq_mode. Options are ILMN_BCL2FASTQ and BCL_PROCESSOR, provided: {}".
                         format(args.fastq_mode))

    # FASTQ input (sample_def)
    hostname = socket.gethostname()
    for idx, sample_def in enumerate(args.sample_def):
        read_path = sample_def.get("read_path")
        if not read_path:
            martian.exit("Must specify a read_path containing FASTQs in each entry of 'sample_def' argument")
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, cellranger-atac does not have permission to open FASTQ folder: %s" % (
                         hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        library_id = sample_def.get("library_id")
        if library_id is not None:
            if not re.match("^[\w-]+$", library_id):
                martian.exit(
                    "Library name may only contain letters, numbers, underscores, and dashes: " + library_id)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not tk_preflight.is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

        if args.fastq_mode == "BCL_PROCESSOR":
            sample_indices, msg = tk_preflight.check_sample_indices(sample_def)
            if sample_indices is None:
                martian.exit(msg)

            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            reads = []
            for sample_index in sample_indices:
                # process interleaved reads
                reads.extend(find_func(read_path, "RA", sample_index, lanes))
            if len(reads) == 0:
                martian.exit("No input FASTQs were found for the requested parameters.")
        elif args.fastq_mode == "ILMN_BCL2FASTQ":
            sample_names = sample_def.get("sample_names", None)
            if sample_names is None:
                martian.exit("Entry {} in sample_def missing required field: sample_names".format(idx))
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            reads1 = []
            reads2 = []
            for sample_name in sample_names:
                r1 = find_func(read_path, BCL2FASTQ_SEQNAMES["R1"], sample_name, lanes)
                r2 = find_func(read_path, BCL2FASTQ_SEQNAMES["R2"], sample_name, lanes)
                if len(r1) != len(r2):
                    martian.exit("Entry {} in sample_defs are missing input FASTQs.".format(idx))
                reads1.extend(r1)
                reads2.extend(r2)
            if len(reads1) == 0 and len(reads2) == 0:
                martian.exit("No input FASTQs were found for the requested parameters.")
        else:
            martian.exit("Unrecognized fastq_mode: {}".format(args.fastq_mode))

    # trim_def['R1'] and ['R2'] must be identical.
    if args.trim_def is not None:
        if len(args.trim_def) == 0:
            martian.exit("trim_def must be a non-empty dictionary.")
        if "R1" not in args.trim_def or "R2" not in args.trim_def:
            martian.exit("trim_def must have R1, R2 fields.")
        if args.trim_def["R1"] != args.trim_def["R2"]:
            martian.exit("trim_def['R1'] and trim_def['R2'] must be identical.")

    # factorization.
    check_factorization(args.factorization)

    # # Reference
    # ref directory structure and timestamps
    ok, msg = check_refdata(args.reference_path, max_contigs=None)
    if ok:
        martian.log_info(msg)
    else:
        martian.exit(msg)

    # usability and format check
    check_reference_format(args.reference_path)

    # Open file handles limit
    if args.check_executables:
        check_filehandle_limit()

    martian.log_info(tk_preflight.record_package_versions())
