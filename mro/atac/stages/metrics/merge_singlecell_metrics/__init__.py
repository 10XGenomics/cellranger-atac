"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Pools single-cell analyses from other stages and calculates bulk metrics from the pooled dataset.
"""
from __future__ import division

import json
import pandas as pd
import numpy as np

from tenkit.stats import robust_divide
from tools import ReferenceManager
from constants import NO_BARCODE
import martian
from metrics.generate_metrics import add_doublet_rate_metrics, add_purity_metrics, add_bulk_targeting_metrics, add_singlecell_sensitivity_metrics


__MRO__ = """
stage MERGE_SINGLECELL_METRICS(
    in  string reference_path,
    in  csv    singlecell_mapping,
    in  csv    singlecell_targets,
    in  csv    singlecell_cells,
    out csv    singlecell,
    out json   summary,
    src py     "stages/metrics/merge_singlecell_metrics",
) using (
    mem_gb = 8,
)
"""

def split(args):
    raise Exception("No split")


def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("No join")

def add_bulk_mapping_metrics(summary_info, singlecell_df, species_list):
    """Waste metric calculations"""

    WASTE_KEYS = [
        "no_barcode",
        "non_cell_barcode",
        "mitochondrial",
        "unmapped",
        "lowmapq",
        "chimeric",
        "duplicate",
        "total"]

    waste_totals = {}
    cell_mask = singlecell_df['cell_id'] != 'None'
    valid_noncell_mask = (singlecell_df['cell_id'] == 'None') & (singlecell_df['barcode'] != NO_BARCODE)

    for key in WASTE_KEYS:
        if key == "no_barcode":
            waste_totals[key] = 0 if all(singlecell_df['barcode'] != 'NO_BARCODE') else\
                singlecell_df['total'].loc[singlecell_df['barcode'] == NO_BARCODE].values[0]
        elif key == "non_cell_barcode":
            waste_totals[key] = singlecell_df['total'][valid_noncell_mask].sum()
        elif key == "total":
            waste_totals[key] = singlecell_df['total'].sum()
        else:
            waste_totals[key] = float(singlecell_df[cell_mask][key].sum())
    for key in WASTE_KEYS:
        if key == "total":
            continue
        summary_info["waste_%s_fragments" % key] = waste_totals[key]
        summary_info["frac_waste_%s" % key] = robust_divide(waste_totals[key], waste_totals["total"])

    if waste_totals['total'] == 0:
        martian.exit("No read pairs were detected in the library. This is likely because \
                     of incorrect input sequenced data potentially due to bioinformatic \
                     mixup. Further execution will be halted.")

    # Special case for frac_waste_no_barcode.
    # We will use its "positive wording" counterpart in the websummary.
    summary_info["frac_valid_barcode"] = 1 - summary_info["frac_waste_no_barcode"]

    cell_mito = singlecell_df[cell_mask]['mitochondrial'].sum()
    cell_unmapped = singlecell_df[cell_mask]['unmapped'].sum()
    cell_lowmapq = singlecell_df[cell_mask]['lowmapq'].sum()
    cell_dup = singlecell_df[cell_mask]['duplicate'].sum()
    cell_passed = singlecell_df[cell_mask]['passed_filters'].sum()

    summary_info['waste_cell_mito_fragments'] = cell_mito
    summary_info['waste_cell_unmapped_fragments'] = cell_unmapped
    summary_info['waste_cell_lowmapq_fragments'] = cell_lowmapq
    summary_info['waste_cell_dup_fragments'] = cell_dup
    summary_info['total_usable_fragments'] = cell_passed

    summary_info['waste_ratio_mito_cells'] = robust_divide(cell_mito, cell_dup + cell_passed)
    summary_info['waste_ratio_unmapped_cells'] = robust_divide(cell_unmapped, cell_dup + cell_passed)
    summary_info['waste_ratio_lowmapq_cells'] = robust_divide(cell_lowmapq, cell_dup + cell_passed)
    summary_info['waste_ratio_dup_cells'] = robust_divide(cell_dup, cell_passed)

    noncell_mito = singlecell_df[valid_noncell_mask]['mitochondrial'].sum()
    noncell_unmapped = singlecell_df[valid_noncell_mask]['unmapped'].sum()
    noncell_lowmapq = singlecell_df[valid_noncell_mask]['lowmapq'].sum()
    noncell_dup = singlecell_df[valid_noncell_mask]['duplicate'].sum()
    noncell_passed = singlecell_df[valid_noncell_mask]['passed_filters'].sum()
    summary_info['total_passed_filter_fragments'] = cell_passed + noncell_passed

    summary_info['waste_ratio_mito_noncells'] = robust_divide(noncell_mito, noncell_dup + noncell_passed)
    summary_info['waste_ratio_unmapped_noncells'] = robust_divide(noncell_unmapped, noncell_dup + noncell_passed)
    summary_info['waste_ratio_lowmapq_noncells'] = robust_divide(noncell_lowmapq, noncell_dup + noncell_passed)
    summary_info['waste_ratio_dup_noncells'] = robust_divide(noncell_dup, noncell_passed)

    cell_barcodes_total = singlecell_df['total'][cell_mask].sum()
    noncell_barcodes_total = singlecell_df['total'][valid_noncell_mask].sum()
    summary_info['num_fragments'] = waste_totals['total']
    summary_info['frac_waste_overall_nondup'] = robust_divide(summary_info['num_fragments'] - (cell_dup + cell_passed), summary_info['num_fragments'])
    summary_info['frac_waste_cell_nondup'] = robust_divide(cell_barcodes_total - (cell_dup + cell_passed), cell_barcodes_total)
    summary_info['frac_valid_noncell'] = robust_divide(noncell_barcodes_total, noncell_barcodes_total + cell_barcodes_total)
    # Next metric is similar to "Reads Mapped Confidently to Genome" in cellranger.
    mapped_confidently = noncell_dup + noncell_passed + cell_dup + cell_passed
    summary_info['frac_mapped_confidently'] = robust_divide(mapped_confidently, summary_info['num_fragments'])

    if not np.isnan(summary_info['frac_mapped_confidently']):
        unmapped_frac = 1 - summary_info['frac_mapped_confidently']
        if unmapped_frac > 0.90:
            martian.exit("%.1f %% of read pairs were not mapped to the supplied \
reference genome. This is likely the consequence of a sample \
mixup or very low sequencing quality. Further execution will \
be halted." % (unmapped_frac * 100))

    summary_info['waste_total_fragments'] = summary_info['num_fragments'] - cell_passed
    summary_info['frac_waste_total'] = robust_divide(summary_info["waste_total_fragments"], summary_info['num_fragments'])

    return summary_info

def main(args, outs):
    if args.singlecell_mapping is None or args.singlecell_targets is None or args.singlecell_cells is None:
        outs.singlecell = None
        outs.summary = None
        return

    ref = ReferenceManager(args.reference_path)
    species_list = ref.list_species()

    # Merge the input singlecell data into a single dataframe and write it out
    mapping = pd.read_csv(args.singlecell_mapping)
    cells = pd.read_csv(args.singlecell_cells)
    targeting = pd.read_csv(args.singlecell_targets)

    merged = mapping.merge(cells, how="left", on="barcode", sort=False, validate="one_to_one")
    merged["cell_id"] = merged["cell_id"].fillna("None")
    for column in merged.columns:
        if column.endswith("_cell_barcode") or column.startswith("passed_filters_") or column.startswith("peak_region_fragments_"):
            merged[column] = merged[column].fillna(0).astype(int)

    merged = merged.merge(targeting, how="left", on="barcode", sort=False, validate="one_to_one")
    keys = ["{}_fragments".format(region)
            for region in ["TSS", "DNase_sensitive_region", "enhancer_region",
                           "promoter_region", "on_target", "blacklist_region", "peak_region"]] + ["peak_region_cutsites"]
    for column in keys:
        merged[column] = merged[column].fillna(0).astype(int)
    merged.to_csv(outs.singlecell, index=None)

    summary_info = {}

    summary_info = add_bulk_targeting_metrics(summary_info, merged, species_list)
    summary_info = add_doublet_rate_metrics(summary_info, merged, species_list)
    summary_info = add_purity_metrics(summary_info, merged, species_list)
    summary_info = add_bulk_mapping_metrics(summary_info, merged, species_list)
    summary_info = add_singlecell_sensitivity_metrics(summary_info, merged, species_list)

    with open(outs.summary, 'w') as summary_file:
        summary_file.write(json.dumps(summary_info, indent=4))
