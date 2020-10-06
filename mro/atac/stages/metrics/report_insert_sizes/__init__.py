"""
Generate a CSV of insert size distributions for each barcode.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division

import json
from collections import Counter, OrderedDict

import martian
import numpy as np
from scipy import signal, stats

import utils
from tools import ReferenceManager, open_fragment_file, combine_csv

__MRO__ = """
stage REPORT_INSERT_SIZES(
    in  tsv.gz fragments,
    in  bool   exclude_non_nuclear,
    in  string reference_path,
    out csv    insert_sizes,
    out json   insert_summary,
    src py     "stages/metrics/report_insert_sizes",
) split (
    in  file   barcode,
    out file   total,
)
"""

MAX_INSERT_SIZE = 1000
GT_MAX_INSERT_SIZE = '>{}'.format(MAX_INSERT_SIZE)


def gen_yvals(xvals, counts, norm_linear=False):
    yvals = np.array([counts[x - 1] for x in xvals])
    if not norm_linear:
        return yvals
    slope, intercept, _, _, _ = stats.linregress(xvals, yvals)
    return yvals - (slope * xvals + intercept)


def find_peak(xvals, counts, target, maxdist):
    """Look for a frequency peak around a given target"""
    try:
        freqs = np.fft.rfftfreq(len(xvals))
        amps = abs(np.fft.rfft(gen_yvals(xvals, counts, True)))
        peak_freqs = freqs[signal.find_peaks_cwt(amps, np.arange(1, 25))]
        diffs = abs((1 / peak_freqs) - target)
        idx = np.argmin(diffs)
        return 1 / peak_freqs[idx] if diffs[idx] <= maxdist else None
    except IndexError:
        # In the case that counts are zero for given xvals
        return None


def split(args):
    if args.fragments is None:
        return {'chunks': [], 'join': {}}

    # as the fragments file is not sorted by barcodes, we iterate through the files to get a list of ordered bcs
    barcodes = list({bc for _, _, _, bc, _ in open_fragment_file(args.fragments)})

    # chunk on barcodes
    barcode_chunks = utils.get_chunks(len(barcodes), 30)
    chunks = []
    for num, bc_chunk in enumerate(barcode_chunks):
        bc_path = martian.make_path('barcode_{}.txt'.format(num))
        with open(bc_path, 'w') as f:
            f.write('\n'.join(barcodes[bc_chunk[0]: bc_chunk[1]]))
        chunks.append({'barcodes': bc_path})
    return {'chunks': chunks, 'join': {'__mem_gb': 16}}

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    with open(args.barcodes, 'r') as barcode_file:
        barcodes_dict = OrderedDict((bc.strip('\n'), num) for num, bc in enumerate(barcode_file))
    outs.insert_summary = None

    if args.fragments is None or len(barcodes_dict) == 0:
        outs.insert_sizes = None
        outs.total = None
        return

    ref_contig_manager = ReferenceManager(args.reference_path)

    # iterate over fragments and count fragment sizes for each barcode
    insert_sizes = {bc: Counter() for bc in barcodes_dict.iterkeys()}
    primary_contigs = set(ref_contig_manager.primary_contigs(allow_sex_chromosomes=True))
    for contig, start, stop, barcode, _ in open_fragment_file(args.fragments):
        if barcode not in barcodes_dict:
            continue
        if args.exclude_non_nuclear and contig not in primary_contigs:
            continue
        size = stop - start
        insert_sizes[barcode][str(size) if size <= MAX_INSERT_SIZE else GT_MAX_INSERT_SIZE] += 1

    # compute total and write out csv
    total = np.zeros(MAX_INSERT_SIZE)
    with open(outs.insert_sizes, 'w') as outfile:
        outfile.write(','.join(
            ['Barcode'] + [str(n) for n in range(1, MAX_INSERT_SIZE + 1)] + ['>{}'.format(MAX_INSERT_SIZE)]) + '\n')
        for barcode in insert_sizes:
            outfile.write(','.join([barcode] +
                                   [str(insert_sizes[barcode][str(n)]) for n in range(1, MAX_INSERT_SIZE + 1)] +
                                   [str(insert_sizes[barcode][GT_MAX_INSERT_SIZE])]) + '\n')
            for n in range(1, 1001):
                total[n - 1] += insert_sizes[barcode][str(n)]

    # write out totals for reduce in join
    np.savetxt(outs.total, total, delimiter=',')

def join(args, outs, chunk_defs, chunk_outs):
    if args.fragments is None:
        outs.insert_sizes = None
        outs.insert_summary = None
        return

    # combine insert_sizes
    combine_csv([chunk_out.insert_sizes for chunk_out in chunk_outs if chunk_out.insert_sizes is not None],
                outs.insert_sizes, header_lines=1)

    # compute total
    total = np.zeros(MAX_INSERT_SIZE)
    for chunk_out in chunk_outs:
        if chunk_out.total is not None:
            total += np.genfromtxt(chunk_out.total, 'float64')

    # NOTE: this assumes MAX_INSERT_SIZE > 701
    xvals = {'full': np.arange(1, MAX_INSERT_SIZE + 1),
             'twist': np.arange(50, 401),
             'nucleosome': np.arange(150, 701)}

    # Calculate metrics
    twist_period = find_peak(xvals['twist'], total, 12, 5)
    nucleosome_period = find_peak(xvals['nucleosome'], total, 180, 30)

    yvals = gen_yvals(xvals['full'], total)
    non_nucleosomal = yvals[xvals['full'] < 147].sum()
    nuc1 = yvals[(147 <= xvals['full']) & (xvals['full'] <= 147 * 2)].sum()

    summary_data = {
        'insert_twist_period': twist_period,
        'insert_nucleosome_period': nucleosome_period,
        'frac_fragments_nfr': non_nucleosomal / yvals.sum(),
        'frac_fragments_nuc': nuc1 / yvals.sum(),
        'frac_fragments_nfr_or_nuc': (non_nucleosomal + nuc1) / yvals.sum(),
    }
    with open(outs.insert_summary, 'w') as outfile:
        outfile.write(json.dumps(summary_data, indent=4))
