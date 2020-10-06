#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#

import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import martian
import subprocess
import json
import os
import numpy as np
from pybedtools import BedTool
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_json
from tools.ref_manager import ReferenceManager
from utils import generate_genome_tag
from constants import TRANSCRIPT_ANNOTATION_GENE_TYPES

MIN_MEM_GB = 4

__MRO__ = """
stage CLOUPE_PREPROCESS(
    in  string     pipestance_type,
    in  string     sample_id,
    in  string     sample_desc,
    in  string     reference_path,
    in  h5         analysis,
    in  h5         feature_barcode_matrix,
    in  bed        peaks,
    in  tsv.gz.tbi fragments_index,
    in  json       metrics_json,
    in  csv        aggregation_csv,
    in  json       gem_group_index_json,
    in  bool       no_secondary_analysis,
    out cloupe     output_for_cloupe,
    out json       gem_group_index_json,
    src py         "stages/cloupe/atac_cloupe_preprocess",
) split (
)
"""

def get_gem_group_index_json(args, outs):
    if args.gem_group_index_json:
        cr_io.copy(args.gem_group_index_json, outs.gem_group_index_json)
    else:
        generated_index = cr_matrix.get_gem_group_index(args.feature_barcode_matrix)
        if generated_index is not None:
            with open(outs.gem_group_index_json, 'w') as outfile:
                tk_json.dump_numpy({"gem_group_index": generated_index}, outfile)
        else:
            outs.gem_group_index_json = None
    return outs.gem_group_index_json

def get_annotation_gene_types(args):
    """
    Return the gene types to use to filter genes/transcript
    annotations by.
    """
    ref_mgr = ReferenceManager(args.reference_path)
    tss = BedTool(ref_mgr.tss_track)
    if tss.field_count() == 7:
        return TRANSCRIPT_ANNOTATION_GENE_TYPES
    else:
        return None

def do_not_make_cloupe(args):
    """
    Returns True if there is a reason why this stage should not attempt to
    generate a .cloupe file
    """
    if args.no_secondary_analysis:
        martian.log_info("Skipping .cloupe generation by instruction (--no-secondary-analysis)")
        return True
    if args.analysis is None or not os.path.exists(args.analysis):
        martian.log_info("Skipping .cloupe generation due to missing analysis hdf5 file")
        return True
    if args.feature_barcode_matrix is None or not os.path.exists(args.feature_barcode_matrix):
        martian.log_info("Skipping .cloupe generation due to missing or zero-length feature-barcode matrix")
        return True
    ref_mgr = ReferenceManager(args.reference_path)
    if len(ref_mgr.list_species()) > 1:
        martian.log_info("Skipping .cloupe generation as the sample is composed of multiple genomes")
        return True
    return False

def split(args):
    # no mem usage if skipped
    if do_not_make_cloupe(args):
        return {'chunks': [], 'join': {'__mem_gb': MIN_MEM_GB}}

    # Multiplicative factor:
    # 2 for storing csc + csr matrix, 0.5 for conversion, 0.4 for the dense propZ matrix, 0.5 for safety
    matrix_mem_gb = int(np.ceil(8 + 3.4 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.feature_barcode_matrix)))
    join = {'__mem_gb': max(matrix_mem_gb, MIN_MEM_GB)}
    return {'chunks': [], 'join': join}


def main(args, outs):
    martian.throw("No chunks defined")


def get_contig_info(reference_path):
    manager = ReferenceManager(reference_path)
    contig_info = {
        "contig_order": {},
        "contig_lengths": {}
    }
    contig_lengths = manager.get_contig_lengths()
    for idx, contig in enumerate(manager.contigs["primary_contigs"]):
        contig_info["contig_order"][contig] = idx
        contig_info["contig_lengths"][contig] = contig_lengths[contig]

    genomes = generate_genome_tag(reference_path)
    contig_info["species"] = genomes
    return contig_info


def join(args, outs, chunk_defs, chunk_outs):
    if do_not_make_cloupe(args):
        outs.output_for_cloupe = None
        return

    reference = ReferenceManager(args.reference_path)

    contig_info_fn = martian.make_path("contig_info.json")
    with open(contig_info_fn, 'w') as outfile:
        contig_info = get_contig_info(args.reference_path)
        json.dump(contig_info, outfile)

    gem_group_index_json = get_gem_group_index_json(args, outs)

    call = ["crconverter",
            args.sample_id,
            args.pipestance_type,
            "--matrix", args.feature_barcode_matrix,
            "--analysis", args.analysis,
            "--output", outs.output_for_cloupe,
            "--description", '"' + args.sample_desc + '"',
            "--peaks", args.peaks,
            "--fragmentsindex", args.fragments_index,
            "--geneannotations", reference.genes,
            "--contiginfo", contig_info_fn,
            ]

    if args.metrics_json is not None:
        call.extend(["--metrics", args.metrics_json])
    if args.aggregation_csv is not None:
        call.extend(["--aggregation", args.aggregation_csv])
    if gem_group_index_json is not None:
        call.extend(["--gemgroups", gem_group_index_json])
    transcript_gene_types = get_annotation_gene_types(args)
    if transcript_gene_types is not None:
        call.extend(["--geneannotationtypes", ",".join(transcript_gene_types)])

    # the sample desc may be unicode, so send the whole
    # set of args str utf-8 to check_output
    unicode_call = [arg.encode('utf-8') for arg in call]

    # but keep the arg 'call' here because log_info inherently
    # attempts to encode the message... (TODO: should log_info
    # figure out the encoding of the input string)
    martian.log_info("Running crconverter: %s" % " ".join(call))
    try:
        results = tk_subproc.check_output(unicode_call)
        martian.log_info("crconverter output: %s" % results)
    except subprocess.CalledProcessError as e:
        outs.output_for_cloupe = None
        martian.throw("Could not generate .cloupe file: \n%s" % e.output)
