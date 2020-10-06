
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

aggregation metrics are devised from final fragments file
"""

import martian
from tools.io import get_counts_by_barcode
from tools import ReferenceManager
from metrics.generate_metrics import (add_doublet_rate_metrics, add_purity_metrics, add_bulk_targeting_metrics,
                                      add_singlecell_sensitivity_metrics)
from utils import get_cell_barcodes
import numpy as np

import json
import pickle
import pandas as pd

__MRO__ = """
stage REGENERATE_ANALYSIS_METRICS(
    in  string     reference_path,
    in  bed        peaks,
    in  csv        cell_barcodes,
    in  tsv.gz     fragments,
    in  tsv.gz.tbi fragments_index,
    out json       regenerated_metrics,
    out csv        singlecell,
    out csv        cell_barcodes,
    src py         "stages/common/regenerate_analysis_metrics",
) split (
    in  string     contig,
    out pickle     target_counts_by_barcode,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


def split(args):
    if args.fragments is None:
        return {'chunks': [], 'join': {}}

    if args.peaks is None:
        martian.throw("peaks BED file expected")
    if args.cell_barcodes is None:
        martian.throw("cell barcodes CSV file expected")

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    chunks = []
    for contig in all_contigs:
        chunks.append({'contig': contig,
                       '__mem_gb': 4})

    return {'chunks': chunks, 'join': {'__mem_gb': 8}}


def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    if args.fragments is None:
        outs.regenerated_metrics = None
        outs.singlecell = None
        return

    target_counts_by_barcode = {}
    ref_mgr = ReferenceManager(args.reference_path)
    for chunk_in, chunk_out in zip(chunk_defs, chunk_outs):
        with open(chunk_out.target_counts_by_barcode, 'r') as infile:
            chunk_counts = pickle.load(infile)
        for barcode, barcode_counts in chunk_counts.iteritems():
            if barcode not in target_counts_by_barcode:
                target_counts_by_barcode[barcode] = barcode_counts
            else:
                for key, value in barcode_counts.iteritems():
                    if key == 'cell_id':
                        target_counts_by_barcode[barcode][key] = value
                    else:
                        target_counts_by_barcode[barcode][key] += value

    species_list = ref_mgr.list_species()
    keys = ["{region}_fragments".format(region=reg)
            for reg in ["TSS", "DNase_sensitive_region", "enhancer_region",
                        "promoter_region", "on_target", "blacklist_region", "peak_region"]] +\
           ["peak_region_cutsites", "passed_filters", "duplicate", "cell_id"] +\
           ["is_{}_cell_barcode".format(species) for species in species_list]
    if len(species_list) > 1:
        keys += ["passed_filters_{}".format(species) for species in species_list] +\
                ["peak_region_fragments_{}".format(species) for species in species_list]

    with open(outs.singlecell, 'w') as outfile:
        outfile.write("barcode,")
        outfile.write(",".join(keys))
        outfile.write("\n")
        for barcode in sorted(target_counts_by_barcode.keys()):
            outfile.write("{},".format(barcode))
            outfile.write(",".join([str(target_counts_by_barcode[barcode][key]) for key in keys]))
            outfile.write("\n")

    # write cell barcodes if uniques > 0 (i.e. subsampling didn't lose barcodes)
    # overwrite the singlecell.csv with update cell calls
    scdf = pd.read_csv(outs.singlecell, sep=',')
    scdf['cell_id'] = np.full(len(scdf), "None")
    ctg_mgr = ReferenceManager(args.reference_path)
    for species in ctg_mgr.list_species():
        species_cell_mask = (scdf['is_{}_cell_barcode'.format(species)] >= 1) & (scdf['passed_filters'] > 0)
        scdf['is_{}_cell_barcode'.format(species)] = np.where(species_cell_mask, 1, 0)
        scdf['cell_id'][species_cell_mask] = np.array(["{}_cell_{}".format(species, num) for num in xrange(np.sum(species_cell_mask))])
    scdf.to_csv(outs.singlecell, sep=',', index=False)
    cell_barcodes = get_cell_barcodes(outs.singlecell, args.reference_path, with_species=True)
    with open(outs.cell_barcodes, 'w') as f:
        for species in cell_barcodes:
            f.write(species + "," + ",".join(cell_barcodes[species]) + "\n")

    # write frag metrics
    summary_info = {}
    summary_info = add_bulk_targeting_metrics(summary_info, scdf, species_list)
    summary_info = add_doublet_rate_metrics(summary_info, scdf, species_list)
    summary_info = add_purity_metrics(summary_info, scdf, species_list)
    summary_info = add_singlecell_sensitivity_metrics(summary_info, scdf, species_list)
    for species in species_list:
        key_suffix = "" if len(species_list) == 1 else "_{}".format(species)
        summary_info['annotated_cells{}'.format(key_suffix)] = scdf['is_{}_cell_barcode'.format(species)].sum()

    with open(outs.regenerated_metrics, 'w') as summary_file:
        summary_file.write(json.dumps(summary_info, indent=4))


def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    if args.fragments is None:
        return

    # create data from frags x barcodes
    _, target_counts, _, _ = get_counts_by_barcode(args.reference_path, args.peaks, args.fragments,
                                                   fragments_index=args.fragments_index, contig=args.contig,
                                                   known_cells=args.cell_barcodes)
    with open(outs.target_counts_by_barcode, 'w') as outfile:
        pickle.dump(target_counts, outfile)
