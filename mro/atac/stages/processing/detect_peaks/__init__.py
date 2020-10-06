"""
Detect open chromatin regions in the genome as measured from ATAC-seq fragment ends.

Models the # of cut sites in a window around each position with a zero-inflated mixture model with:
    - a geometric distribution to model the zero inflation
    - a negative binomial distribution for the background
    - a geometric distribution distribution for signal (peaks)

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
from __future__ import division

import json
import pickle

from tools.regions import sort_and_uniq_bed
from tools import ReferenceManager, combine_csv
from utils import count_bases_in_peaks, cut_site_counter_from_bedgraph
from collections import Counter
from analysis.peaks import MixedModelParams, estimate_final_threshold
from constants import PEAK_MERGE_DISTANCE, PEAK_ODDS_RATIO

__MRO__ = """
stage DETECT_PEAKS(
    in  bedgraph   cut_sites,
    in  path       reference_path,
    in  pickle     count_dict,
    out bed        peaks,
    out json       peak_metrics,
    src py         "stages/processing/detect_peaks",
) split (
    in  string     contig,
    in  float[]    params,
    in  float      threshold,
) using (
    mem_gb = 6,
)
"""

def identify_peaks(site_iter, out_peaks, threshold):
    """Write out all peaks found in the count data to an output file, using a single fixed threshold of windowed
    counts."""
    summary_metrics = Counter()
    summary_metrics['total_peaks_detected'] = 0

    window = None
    with open(out_peaks, 'w') as outfile_peaks:
        for chrom, pos, windowed_count in site_iter:
            if windowed_count < threshold:
                continue
            if window is None:
                window = [chrom, pos, pos + 1]
            elif pos - PEAK_MERGE_DISTANCE <= window[2] and chrom == window[0]:
                window[2] = pos + 1
            else:
                summary_metrics['total_peaks_detected'] += 1
                outfile_peaks.write('{}\t{}\t{}\n'.format(*window))
                window = [chrom, pos, pos + 1]
        if window is not None:
            summary_metrics['total_peaks_detected'] += 1
            outfile_peaks.write('{}\t{}\t{}\n'.format(*window))
    return summary_metrics

def split(args):
    if args.cut_sites is None:
        return {'chunks': [], 'join': {}}

    ctg_mgr = ReferenceManager(args.reference_path)
    all_contigs = ctg_mgr.primary_contigs(allow_sex_chromosomes=True)

    # First we estimate the peak threshold and output a bedgraph of cut site positions and counts
    with open(args.count_dict, 'r') as f:
        count_dict = pickle.load(f)
    threshold, params = estimate_final_threshold(count_dict, PEAK_ODDS_RATIO)

    chunks = []
    for contig in all_contigs:
        chunks.append({'params': params,
                       'contig': contig,
                       'threshold': threshold,
                       '__mem_gb': 2})

    return {'chunks': chunks, 'join': {'__mem_gb': 5}}

def main(args, outs):
    '''Call peaks on chromosomes'''
    if args.cut_sites is None:
        outs.peak_metrics = None
        outs.peaks = None
        return

    # identify peaks based on the estimated thresholds
    site_iter = cut_site_counter_from_bedgraph(args.cut_sites, threshold=args.threshold, contig=args.contig)
    summary_metrics = identify_peaks(site_iter, outs.peaks, args.threshold)

    with open(outs.peak_metrics, 'w') as outfile:
        json.dump(summary_metrics, outfile)

def join(args, outs, chunk_defs, chunk_outs):
    '''Join peaks'''

    if args.cut_sites is None:
        outs.peak_metrics = None
        outs.peaks = None
        return

    # merge peaks in contig order in chunks
    chunk_peaks = [chunk_out.peaks for chunk_out in chunk_outs if chunk_out.peaks is not None]
    combine_csv(chunk_peaks, outs.peaks, header_lines=0)
    sort_and_uniq_bed(outs.peaks, args.reference_path)

    # Merge metrics by reduction
    peak_metrics = Counter()
    for chunk_out in chunk_outs:
        with open(chunk_out.peak_metrics, 'r') as f:
            peak_metrics += Counter(json.load(f))

    # count bases in peaks
    peak_metrics['bases_in_peaks'] = count_bases_in_peaks(args.reference_path, outs.peaks)

    # add fit data to summary
    fit_threshold = chunk_defs[0].threshold
    if chunk_defs[0].params is not None:
        fit_params = MixedModelParams(*chunk_defs[0].params)
        for name in fit_params._fields:
            peak_metrics["peakcaller_{}".format(name)] = getattr(fit_params, name)
    peak_metrics["peakcaller_threshold"] = fit_threshold
    with open(outs.peak_metrics, 'w') as f:
        json.dump(peak_metrics, f)
