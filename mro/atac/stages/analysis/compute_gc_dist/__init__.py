
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

compute gc distribution in peaks
"""

# Imports
from __future__ import division
from tools import ReferenceManager
from tools.io import peak_reader
import utils
import martian
from constants import VALID_BASES as JASPAR_MOTIF_BASE_ORDER

import pyfasta
from collections import Counter
import pickle
import numpy as np


# MRO docstring
__MRO__ = """
stage COMPUTE_GC_DISTRIBUTION(
    in  bed    peaks,
    in  string reference_path,
    out pickle GCdict,
    src py     "stages/analysis/compute_gc_dist",
) split (
) using (
    volatile = strict,
)
"""

LOW_GC = 0.33
HIGH_GC = 0.7
NBINS = 25

def get_GCbinned_peaks_and_bg(peaks, genome_fa, GCbins, pseudocount=1.0):
    '''Calculate the G-C and A-T % within the peaks for each bin'''

    def findbin(peakGC, GCbins):
        '''Expected sorted GCbins'''
        for gc in GCbins:
            if peakGC < gc[1] and peakGC >= gc[0]:
                return gc
        if abs(peakGC - GCbins[-1][1]) < 1e-5:
            return GCbins[-1]

    GCdict = {}
    base_counter = {}
    for gc in GCbins:
        GCdict[gc] = {}
        GCdict[gc]['counter'] = Counter({base: 0 for base in JASPAR_MOTIF_BASE_ORDER})
        GCdict[gc]['peaks'] = []

    for num, peak in enumerate(peaks):
        peakGC, base_counter = utils.get_peak_GC_counts(peak, genome_fa)
        gc = findbin(peakGC, GCbins)
        GCdict[gc]['peaks'].append(num)
        GCdict[gc]['counter'] += Counter(base_counter)

    for gc in GCbins:
        GCdict[gc]['total_bases'] = sum(GCdict[gc]['counter'].values())

        freq = {}
        for base in JASPAR_MOTIF_BASE_ORDER:
            freq[base] = np.divide(GCdict[gc]['counter'][base] + pseudocount,
                                   GCdict[gc]['total_bases'] + 4 * pseudocount, dtype=float)

        # the freq is computed using only the "+" strand. To get the freq of both strand, average A&T and G&C
        bg = {}
        bg['A'] = (freq['A'] + freq['T']) / 2
        bg['C'] = (freq['G'] + freq['C']) / 2
        bg['G'] = bg['C']
        bg['T'] = bg['A']
        GCdict[gc]['counter'] = [bg[base] for base in JASPAR_MOTIF_BASE_ORDER]

    return GCdict

def split(args):
    ref_mgr = ReferenceManager(args.reference_path)
    return {'chunks': [], 'join': {'__mem_gb': 4, '__vmem_gb': int(np.ceil(ref_mgr.get_vmem_est())) + 3}}

def main(args, outs):
    martian.throw("main not implemented")

def join(args, outs, chunk_defs, chunk_outs):
    """Compute base background in each peak."""
    ref_mgr = ReferenceManager(args.reference_path)
    npeaks = utils.quick_line_count(args.peaks) if args.peaks else 0

    if len(ref_mgr.list_species()) > 1 or npeaks == 0 or ref_mgr.motifs is None:
        outs.GCdist = None
        return

    # get peak-GC distribution
    genome_fa = pyfasta.Fasta(ref_mgr.fasta, key_fn=lambda x: x.split()[0])
    GCdist = [utils.get_peak_GC_counts(peak, genome_fa, counts=False) for peak in peak_reader(args.peaks)]

    # compute base background from peaks in bins
    # merge extreme GC bins with adjoining ones if they're too narrow for motif scanner to work correctly
    GCbounds = []
    nbins = NBINS
    for n, gc in enumerate(np.percentile(GCdist, np.linspace(0, 100, nbins + 1, endpoint=True), interpolation='lower')):
        if n == 0 or n == nbins:
            GCbounds += [gc]
            continue
        if gc >= LOW_GC and gc < HIGH_GC:
            GCbounds += [gc]
    GCbins = sorted(list(set(zip(GCbounds, GCbounds[1:]))))  # uniqify
    peaks = peak_reader(args.peaks)
    GCdict = get_GCbinned_peaks_and_bg(peaks, genome_fa, GCbins)

    # dump
    with open(outs.GCdict, 'w') as f:
        pickle.dump(GCdict, f)
