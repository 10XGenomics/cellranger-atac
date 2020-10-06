#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

import martian
import tenkit.preflight as tk_preflight
from preflights import (check_refdata, check_vmem_for_reference, check_factorization,
                        check_reference_format, check_filehandle_limit, check_aggr_csv, check_sample_id)
from constants import ALLOWED_NORMALIZATIONS

__MRO__ = """
stage ATAC_AGGR_PREFLIGHT(
    in  string   sample_id,
    in  string   reference_path,
    in  csv      aggr_csv,
    in  string   normalization,
    in  string[] factorization,
    in  bool     check_executables,
    src py       "stages/preflight/atac_aggr",
) split (
)
"""

def split(args):
    mem_gb = 6  # HACK: reserve a basic chunk when aggr list is long
    vmem_gb = mem_gb + check_vmem_for_reference(args.reference_path)
    return {'chunks': [], 'join': {'__mem_gb': mem_gb, '__vmem_gb': vmem_gb}}

def main(args, outs):
    martian.throw("main not implemented")

def join(args, outs, chunk_defs, chunk_outs):
    # Sample ID / pipestance name
    check_sample_id(args.sample_id)

    # factorization.
    check_factorization(args.factorization)

    # normalization checks
    if args.normalization not in ALLOWED_NORMALIZATIONS:
        martian.exit("Unsupported normalization method provided. Options are {}.".
                     format(", ".join(ALLOWED_NORMALIZATIONS)))

    # # Reference
    # ref directory structure and timestamps
    ok, msg = check_refdata(args.reference_path, max_contigs=None)
    if ok:
        martian.log_info(msg)
    else:
        martian.exit(msg)

    # usability and check file formats
    check_reference_format(args.reference_path)

    # aggr csv checks
    if args.aggr_csv is None:
        martian.exit("aggregation csv must be provided")
    check_aggr_csv(args.aggr_csv, args.reference_path)

    # Open file handles limit
    if args.check_executables:
        check_filehandle_limit()

    martian.log_info(tk_preflight.record_package_versions())
