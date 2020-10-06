'''
 Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

 Scan the peaks for motif matches using a pval threshold.
'''

from __future__ import division

import martian
from analysis.motifs import Motifs
from tools import ReferenceManager
from tools.io import combine_csv, peak_reader
from constants import PWM_MATCH_PVAL_THRESHOLD
import pickle
import utils
import numpy as np

__MRO__ = """
stage SCAN_MOTIFS(
    in  pickle globalGCdict,
    in  bed    peaks,
    in  string reference_path,
    in  float  pwm_threshold,
    out bed    peak_motif_hits,
) split (
    in  file   GCdict,
) using (
    volatile = strict,
)
"""

def split(args):
    """Compute base background in split and use it in each chunk."""

    ref_mgr = ReferenceManager(args.reference_path)
    npeaks = utils.quick_line_count(args.peaks) if args.peaks else 0
    if len(ref_mgr.list_species()) > 1 or npeaks == 0 or ref_mgr.motifs is None:
        chunk_def = [{'skip': True}]
        return {'chunks': chunk_def}

    with open(args.globalGCdict, 'r') as f:
        GCdict = pickle.load(f)

    GCdict_paths = {}
    GCbins = sorted(GCdict.keys())
    for gc in GCbins:
        GCdict_paths[gc] = martian.make_path('GCdict_{}_{}'.format(gc[0], gc[1]))
        with open(GCdict_paths[gc], 'w') as dump:
            pickle.dump(GCdict[gc], dump)

    # write rows of each chunk to a new peak file
    mem_in_gb = 8
    chunk_def = [{'__mem_gb': mem_in_gb,
                  '__vmem_gb': mem_in_gb + int(np.ceil(ref_mgr.get_vmem_est())) + 1,
                  'skip': False,
                  'GCdict': GCdict_paths[chunk]} for chunk in GCbins]
    return {'chunks': chunk_def}

def main(args, outs):
    if args.skip:
        outs.peak_motif_hits = None
        return

    with open(args.GCdict, 'r') as dump:
        GCdict = pickle.load(dump)

    # annotate peak-motif mappings. N.B. the peaks and motifs are 0-indexed
    peaks_iter = peak_reader(args.peaks, select=GCdict['peaks'])
    motifs = Motifs(args.reference_path, bg=GCdict['counter'])
    use_pwm_threshold = args.pwm_threshold if args.pwm_threshold is not None else PWM_MATCH_PVAL_THRESHOLD
    motifs.scan_motif_from_bed(peaks_iter, out_file=outs.peak_motif_hits, out_format="binary-bed",
                               use_genome_bg=False, pvalue=use_pwm_threshold)


def join(args, outs, chunk_defs, chunk_outs):
    if chunk_defs[0].skip:
        outs.peak_motif_hits = None
        return

    # NOTE: we're okay creating an unsorted bed file
    chunk_peak_motif_hits = [chunk.peak_motif_hits for chunk in chunk_outs if chunk.peak_motif_hits is not None]
    combine_csv(chunk_peak_motif_hits, outs.peak_motif_hits, header_lines=0)
