"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Estimates total library complexity -- both bulk and on a single cell basis -- using total and unique fragments.
"""

from __future__ import division

import json
from collections import Counter, defaultdict
from complexity import (estimate_library_complexity_safe, estimate_complexity_dual, mean_unique_from_total_safe,
                        get_unique_and_total_fragments, downsample_counts, interpolate_downsampled_counts)
from tools import combine_csv
import utils
import numpy as np
import tenkit.stats as tk_stats
from constants import DOWNSAMPLED_READS, RPC_30K, RPC_50K, RPC_10K
from tools import open_fragment_file
import martian
import pandas as pd

__MRO__ = """
stage ESTIMATE_LIBRARY_COMPLEXITY(
    in  json   sequencing_summary,
    in  tsv.gz fragments,
    in  csv    cell_barcodes,
    out json   singlecell_complexity,
    out json   bulk_complexity,
    out json   complexity_summary,
    src py     "stages/metrics/estimate_library_complexity",
) split (
    in  file   barcodes,
) using (
    mem_gb = 6,
)
"""

def split(args):
    if args.fragments is None or args.sequencing_summary is None:
        return {'chunks': [], 'join': {}}

    # as the fragments file is not sorted by barcodes, we iterate through the files to get a list of bcs
    barcodes = list({bc for _, _, _, bc, _ in open_fragment_file(args.fragments)})

    # chunk on barcodes
    barcode_chunks = utils.get_chunks(len(barcodes), 20)
    chunks = []
    for num, bc_chunk in enumerate(barcode_chunks):
        bc_path = martian.make_path('barcode_{}.txt'.format(num))
        with open(bc_path, 'w') as f:
            f.write('\n'.join(barcodes[bc_chunk[0]: bc_chunk[1]]))
        chunks.append({'__mem_gb': 4,
                       'barcodes': bc_path})
    return {'chunks': chunks, 'join': {'__mem_gb': 8}}

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    if args.fragments is None or args.sequencing_summary is None:
        outs.singlecell_complexity = None
        outs.complexity_summary = None
        outs.bulk_complexity = None
        return

    singlecell_histogram = defaultdict(Counter)

    # load cell barcodes
    if args.cell_barcodes:
        cell_barcodes = utils.load_cell_barcodes(args.cell_barcodes)

    # get dup counts for each barcode
    with open(args.barcodes, 'r') as f:
        barcodes = {bc.strip('\n') for bc in f}
    for contig, start, stop, barcode, dup_count in open_fragment_file(args.fragments):
        if barcode in barcodes:
            singlecell_histogram[barcode][dup_count] += 1

    # Write out the singlecell complexity series of uniques at fixed downsampling rates
    downsample_rates = np.linspace(1e-5, 0.99, 100)
    total_fragments_per_cell = {}
    with open(outs.singlecell_complexity, 'w') as sc_complexity_file:
        sc_complexity_file.write("Barcode,")
        # NOTE: the extra key added because we append the no-donwsampling case in downsample_counts
        sc_complexity_file.write(",".join([str(rate) for rate in downsample_rates]) + ",1.0\n")
        for barcode in singlecell_histogram:
            if barcode in cell_barcodes:
                _, total_fragments_per_cell[barcode] = get_unique_and_total_fragments(singlecell_histogram[barcode])
                unique_at_ds, _ = downsample_counts(singlecell_histogram[barcode], downsample_rates)
                sc_complexity_file.write(barcode + ",")
                sc_complexity_file.write(",".join(map(str, unique_at_ds.tolist())) + "\n")
    bulk_histogram = sum((singlecell_histogram[barcode] for barcode in singlecell_histogram), Counter())

    # NOTE: repurpose the summary file and bulk complexity files to avoid creating unneeded new outputs
    with open(outs.complexity_summary, 'w') as outfile:
        outfile.write(json.dumps(total_fragments_per_cell, indent=4))
    with open(outs.bulk_complexity, 'w') as outfile:
        outfile.write(json.dumps(bulk_histogram, indent=4))


def join(args, outs, chunk_defs, chunk_outs):
    if args.fragments is None or args.sequencing_summary is None:
        outs.singlecell_complexity = None
        outs.complexity_summary = None
        outs.bulk_complexity = None
        return

    # sequencing info
    with open(args.sequencing_summary, 'r') as sm:
        num_reads = json.load(sm)['num_reads']

    summary_data = {}
    # # Single cell complexity:
    # merge single cell complexity rows
    combine_csv([chunk_out.singlecell_complexity for chunk_out in chunk_outs if chunk_out.singlecell_complexity is not None],
                outs.singlecell_complexity, header_lines=1)

    # merge total_fragments_per_cell info
    total_fragments_per_cell = {}
    for chunk_out in chunk_outs:
        with open(chunk_out.complexity_summary, 'r') as infile:
            total_fragments_per_cell.update(json.load(infile))

    # calculate median unique fragments at various downsamplings
    downsampled_unique_fragments = pd.read_csv(outs.singlecell_complexity, index_col='Barcode')
    downsample_rates = downsampled_unique_fragments.columns.values.astype(float)

    # if no cell barcodes, skip this part
    num_cells = downsampled_unique_fragments.values.shape[0]
    if num_cells == 0:
        outs.singlecell_complexity = None
    else:
        # plot data
        median_total_array = downsample_rates * np.median(total_fragments_per_cell.values())
        mean_total_array = downsample_rates * np.mean(total_fragments_per_cell.values())
        median_unique_array = np.median(downsampled_unique_fragments.values, axis=0)
        # for visualization, total depth will be reported on a per cell basis, to be consistent with cellranger:
        total_depth_per_cell_array = downsample_rates * num_reads / num_cells

        with open(outs.singlecell_complexity, 'w') as outfile:
            outfile.write(json.dumps({'total_depth': total_depth_per_cell_array.tolist(),
                                      'total': mean_total_array.tolist(),
                                      'unique': median_unique_array.tolist()}, indent=4))

        pool1, pool2, frac_in_pool1 = estimate_complexity_dual(median_total_array, median_unique_array)
        median_estimated_complexity = np.nan if np.isnan(frac_in_pool1) else pool1 + pool2

        if not np.isnan(median_estimated_complexity):
            summary_data.update({
                'median_per_cell_total_library_complexity': median_estimated_complexity,
                'median_per_cell_estimated_saturation': tk_stats.robust_divide(median_unique_array[-1], median_estimated_complexity)}
            )

        # per species RRPC and HQ RPC metrics:
        # load species cell barcodes
        cell_barcodes = {}
        if args.cell_barcodes:
            cell_barcodes = utils.load_cell_barcodes(args.cell_barcodes, with_species=True)

        for species in cell_barcodes.keys():
            species_prefix = "_{}".format(species) if len(cell_barcodes) > 1 else ""
            mean_total_array_sp = downsample_rates * np.mean([total_fragments_per_cell[bc] for bc in cell_barcodes[species]])
            downsampled_unique_fragments_sp = pd.DataFrame(downsampled_unique_fragments, index=cell_barcodes[species])
            median_unique_array_sp = np.median(downsampled_unique_fragments_sp.values, axis=0)
            if (mean_total_array_sp[-1] * 2) >= RPC_50K:
                ds_rate = RPC_50K / (mean_total_array_sp[-1] * 2)
                ds_fragment_count = interpolate_downsampled_counts(mean_total_array_sp, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_HQ_RPC{}'.format(RPC_50K, species_prefix): ds_fragment_count})

            if (mean_total_array_sp[-1] * 2) >= RPC_30K:
                ds_rate = RPC_30K / (mean_total_array_sp[-1] * 2)
                ds_fragment_count = interpolate_downsampled_counts(mean_total_array_sp, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_HQ_RPC{}'.format(RPC_30K, species_prefix): ds_fragment_count})

            if (mean_total_array_sp[-1] * 2) >= RPC_10K:
                ds_rate = RPC_10K / (mean_total_array_sp[-1] * 2)
                ds_fragment_count = interpolate_downsampled_counts(mean_total_array_sp, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_HQ_RPC{}'.format(RPC_10K, species_prefix): ds_fragment_count})

            if (num_reads / num_cells) >= RPC_50K:
                ds_rate = RPC_50K / (num_reads / num_cells)
                ds_fragment_count = interpolate_downsampled_counts(total_depth_per_cell_array, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_RRPC{}'.format(RPC_50K, species_prefix): ds_fragment_count})

            if (num_reads / num_cells) >= RPC_30K:
                ds_rate = RPC_30K / (num_reads / num_cells)
                ds_fragment_count = interpolate_downsampled_counts(total_depth_per_cell_array, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_RRPC{}'.format(RPC_30K, species_prefix): ds_fragment_count})

            if (num_reads / num_cells) >= RPC_10K:
                ds_rate = RPC_10K / (num_reads / num_cells)
                ds_fragment_count = interpolate_downsampled_counts(total_depth_per_cell_array, median_unique_array_sp, ds_rate)
                summary_data.update({'median_per_cell_unique_fragments_at_{}_RRPC{}'.format(RPC_10K, species_prefix): ds_fragment_count})

    # # Bulk complexity:
    # compute bulk complexity metrics
    bulk_histogram = Counter()
    for chunk_out in chunk_outs:
        with open(chunk_out.bulk_complexity, 'r') as infile:
            bulk_histogram += Counter(json.load(infile))

    unique, total = get_unique_and_total_fragments(bulk_histogram)
    bulk_complexity = estimate_library_complexity_safe(bulk_histogram, 'compressed')

    if (total * 2) >= DOWNSAMPLED_READS:
        ds_rate = DOWNSAMPLED_READS / (total * 2)
        downsampled_unique, _ = downsample_counts(bulk_histogram, [ds_rate])
        ds_fragment_count = downsampled_unique[0]
        summary_data.update({'bulk_unique_fragments_at_{}_reads'.format(DOWNSAMPLED_READS): ds_fragment_count})

    # if complexity estimation is successful, compute some metrics and fitted plot
    if not np.isnan(bulk_complexity['fraction_pool1']):
        bulk_complexity['plot_total'] = list(np.linspace(1e5, max(bulk_complexity['total']) * 1.2, 100))
        bulk_complexity['plot_estimated_unique'] = list(mean_unique_from_total_safe(bulk_complexity))

        summary_data.update({
            'bulk_total_library_complexity': bulk_complexity['estimated_complexity'],
            'bulk_estimated_saturation': tk_stats.robust_divide(unique, bulk_complexity["estimated_complexity"])
        })

    with open(outs.bulk_complexity, 'w') as outfile:
        outfile.write(json.dumps(bulk_complexity, indent=4))

    if summary_data:
        with open(outs.complexity_summary, 'w') as outfile:
            outfile.write(json.dumps(summary_data, indent=4))
    else:
        outs.complexity_summary = None
