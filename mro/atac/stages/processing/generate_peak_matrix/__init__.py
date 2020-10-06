"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Generate a matrix of barcodes by peaks, with the counts of the # of transposition sites in each
"""

from __future__ import division

import os
from scipy.sparse import csc_matrix, coo_matrix, hstack
from collections import OrderedDict, Counter

import utils
import tenkit.bio_io as tk_bio
import cellranger.matrix as cr_matrix
import cellranger.atac.feature_ref as atac_feature_ref
import cellranger.atac.matrix as atac_matrix
import cellranger.library_constants as cr_lib_constants

from tools import open_fragment_file
import martian

__MRO__ = """
stage GENERATE_PEAK_MATRIX(
    in  string reference_path,
    in  tsv.gz fragments,
    in  bed    peaks,
    out path   raw_matrix_mex,
    out h5     raw_matrix,
    src py     "stages/processing/generate_peak_matrix",
) split (
    in  file   barcodes,
) using (
    mem_gb = 4,
)
"""


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

def join(args, outs, chunk_defs, chunk_outs):
    if args.fragments is None:
        outs.raw_matrix = None
        outs.raw_matrix_mex = None
        return

    # Hstack barcodes to generate full peak matrix
    barcodes = []
    sp_matrix = None
    for i, chunk in enumerate(chunk_outs):
        if chunk.raw_matrix is not None and os.path.exists(chunk.raw_matrix):
            cpm = cr_matrix.CountMatrix.load_h5_file(chunk.raw_matrix)
            if i == 0:
                sp_matrix = cpm.m
            else:
                sp_matrix = hstack([sp_matrix, cpm.m])
            barcodes.extend(cpm.bcs)

    genomes = utils.generate_genome_tag(args.reference_path)
    peaks_def = atac_feature_ref.from_peaks_bed(args.peaks, genomes)
    raw_matrix = cr_matrix.CountMatrix(peaks_def, barcodes, sp_matrix)
    raw_matrix.save_h5_file(outs.raw_matrix, sw_version=martian.get_pipelines_version())
    if not os.path.exists(outs.raw_matrix_mex):
        os.mkdir(outs.raw_matrix_mex)
    atac_matrix.save_mex(raw_matrix, outs.raw_matrix_mex,
                         cr_lib_constants.ATACSEQ_LIBRARY_TYPE,
                         martian.get_pipelines_version())

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    outs.raw_matrix_mex = None
    if args.fragments is None:
        outs.raw_matrix = None
        return

    with open(args.peaks, 'r') as infile:
        full_peaks = tk_bio.get_target_regions(infile)
    with open(args.peaks, 'r') as pfile:
        peaks_dict = OrderedDict(("{}:{}-{}".format(*peak.strip("\n").split("\t")), num) for num, peak in enumerate(pfile))

    with open(args.barcodes, 'r') as barcode_file:
        barcodes_dict = OrderedDict((bc.strip('\n'), num) for num, bc in enumerate(barcode_file))

    if len(barcodes_dict) == 0:
        outs.raw_matrix = None
        return

    # get matrix counts
    peak_bc_counts = Counter()
    for contig, start, stop, barcode, _ in open_fragment_file(args.fragments):
        if barcode not in barcodes_dict:
            continue
        for pos in (start, stop):
            if contig in full_peaks.keys():
                peak = full_peaks[contig].get_region_containing_point(pos)
                if peak is not None:
                    peak_bc_counts[barcodes_dict[barcode], peaks_dict['{}:{}-{}'.format(contig, peak[0], peak[1])]] += 1

    data, col, row = (), (), ()
    if len(peak_bc_counts) > 0:
        data, col, row = zip(*[(val, key[0], key[1]) for key, val in peak_bc_counts.iteritems()])
    sp_matrix = csc_matrix(coo_matrix((data, (row, col)), shape=(len(peaks_dict), len(barcodes_dict)), dtype=int))

    # save as a CountMatrix
    genomes = utils.generate_genome_tag(args.reference_path)
    peaks_def = atac_feature_ref.from_peaks_bed(args.peaks, genomes)
    raw_matrix = cr_matrix.CountMatrix(peaks_def, barcodes_dict.keys(), sp_matrix)
    raw_matrix.save_h5_file(outs.raw_matrix, sw_version=martian.get_pipelines_version())
