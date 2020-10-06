"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
# Generate barcode-TF matrix based on the peak-motif annotation

"""

from __future__ import division

import os
import scipy.sparse as sp
import numpy as np
from collections import OrderedDict, Counter
from statsmodels import robust
from pybedtools import BedTool

import cellranger.atac.feature_ref as atac_feature_ref
import cellranger.atac.matrix as atac_matrix
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.matrix as cr_matrix
import cellranger.library_constants as cr_lib_constants

import martian
from analysis.motifs import Motifs
import utils
from tools import ReferenceManager

__MRO__ = """
stage GENERATE_TF_MATRIX(
    in  path reference_path,
    in  bed  peaks,
    in  bed  peak_motif_hits,
    in  h5   filtered_matrix,
    out h5   filtered_tf_bc_matrix,
    out path filtered_tf_bc_matrix_mex,
    out gz   tf_propZ_matrix,
    src py   "stages/analysis/generate_tf_matrix",
) split (
)
"""


def motifscan_bed_to_sparse_matrix(scan_bed, peak_idx, motif_idx, format='binary'):
    """

    :param scan_bed: BedTool object of the peak_motif_hits.bed from the ANNOTATE_PEAKS output
    :param peak_idx: an OrderedDict of peak coordinates to integer mapping
    :param motif_idx: an OrderedDict of motif name to integer mapping
    :param format: output of the sparse matrix is either binary or count
    :return: a scipy.sparse.csc_matrix
    """
    assert format in ['count', 'binary']

    tf_peak_counter = Counter()
    peak_coor = []
    motif_coor = []
    values = []
    for row in scan_bed:
        peak_id = peak_idx[str('_'.join(row[:3]))]
        motif_id = motif_idx[str(row[3])]

        if format == 'count':
            tf_peak_counter[(peak_id, motif_id)] += 1
        elif format == 'binary':
            tf_peak_counter[(peak_id, motif_id)] = 1

    for key, v in tf_peak_counter.items():
        peak_coor.append(key[0])
        motif_coor.append(key[1])
        values.append(v)

    return peak_coor, motif_coor, values


def _get_peak_indexes(peaks):
    """
    :param peaks: BedTool object
    """
    out_dict = OrderedDict()
    for i, peak in enumerate(peaks):
        peak_id = str('_'.join(peak[:3]))
        out_dict[peak_id] = i
    return out_dict, i + 1


def _get_motif_indexes(motifs):
    """
    :param motif: Motifs object
    """
    out_dict = OrderedDict()
    for i, motif in enumerate(motifs.all_motifs):
        motif_id = motif.name
        out_dict[motif_id] = i
    return out_dict, i + 1

def MADzscore(matrix, axis=1):
    '''Expects a numpy matrix and returns a zscore matrix. This is robust to 1-D outliers as
    it normalizes the distance from the median using the median absolute distance from the median.
    More: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm#Z-Scores'''
    assert matrix.ndim == 2
    medians = np.median(matrix, axis=axis, keepdims=True)
    MAD = np.expand_dims(robust.mad(matrix, axis=axis), axis=axis)
    return (matrix - medians) / MAD

def split(args):
    ref_mgr = ReferenceManager(args.reference_path)
    if args.filtered_matrix is None or args.peak_motif_hits is None or len(ref_mgr.list_species()) > 1:
        return {'chunks': []}

    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix)
    npeaks, nbcs, nnz = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    # assume we will never test more than 1000 TFs and
    # the relative hit-rate of a TF is a generous 1 out of every 10 peaks
    MAX_TF_COUNT = 1000
    MAX_TF_PEAK_SPARSITY = 0.1
    BYTES_PER_INT = np.dtype(int).itemsize
    BYTES_PER_FLOAT = np.dtype(float).itemsize
    BYTES_PER_GB = 1024**3
    ENTRIES_PER_VAL = 3
    predicted_tf_peak_matrix_mem_gb = ENTRIES_PER_VAL * MAX_TF_PEAK_SPARSITY * npeaks * MAX_TF_COUNT * BYTES_PER_INT / BYTES_PER_GB
    predicted_tf_matrix_mem_gb = ENTRIES_PER_VAL * nbcs * MAX_TF_COUNT * BYTES_PER_INT / BYTES_PER_GB
    predicted_tf_propZ_matrix_mem_gb = ENTRIES_PER_VAL * nbcs * MAX_TF_COUNT * BYTES_PER_FLOAT / BYTES_PER_GB
    chunk_mem_gb = int(np.ceil(max(matrix_mem_gb +
                               predicted_tf_peak_matrix_mem_gb * 2 +
                               predicted_tf_matrix_mem_gb * 2 +
                               predicted_tf_propZ_matrix_mem_gb * 2,
                               h5_constants.MIN_MEM_GB)))
    vmem_peak_motif_hits = int(np.ceil(predicted_tf_peak_matrix_mem_gb) * 3 + predicted_tf_peak_matrix_mem_gb)

    # HACK - give big jobs more threads in order to avoid overloading a node
    threads = cr_io.get_thread_request_from_mem_gb(chunk_mem_gb)

    return {'chunks': [],
            'join': {
                '__mem_gb': chunk_mem_gb,
                '__vmem_gb': chunk_mem_gb + vmem_peak_motif_hits + 1,
                '__threads': threads}
            }


def main(args, outs):
    martian.throw('No chunks defined')


def join(args, outs, chunk_defs, chunk_outs):
    ref_mgr = ReferenceManager(args.reference_path)
    if args.filtered_matrix is None or args.peak_motif_hits is None or len(ref_mgr.list_species()) > 1:
        outs.filtered_tf_bc_matrix = None
        outs.filtered_tf_bc_matrix_mex = None
        outs.tf_propZ_matrix = None
        return

    # motif scan is completed in ANNOTATE_PEAKS

    peaks = BedTool(args.peaks)
    motifs = Motifs(args.reference_path)

    peak_motif_hits = BedTool(args.peak_motif_hits)

    # extract peak coordinate to numerical index map
    peak_idx, n_peaks = _get_peak_indexes(peaks)

    # extract motif names to numerical index map
    motif_idx, n_motifs = _get_motif_indexes(motifs)

    # extract 3 lists: peak indexes, motif indexes and counts, each entry correspond to a peak-motif pair
    peak_coor, motif_coor, values = motifscan_bed_to_sparse_matrix(peak_motif_hits, peak_idx, motif_idx, format='binary')

    # convert it to a sparse matrix, default is binary format, motifs are rows and peaks are columns
    tf_peak_matrix = sp.csr_matrix((values, (motif_coor, peak_coor)), shape=(n_motifs, n_peaks), dtype='int32')

    # compute the motif-BC matrix via pooling
    # The current method simply counts the number of hits for a motif inside the peaks in a barcode
    # cast as a CountMatrix
    peak_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrix)
    motif_names = motif_idx.keys()
    barcodes = peak_matrix.bcs
    genomes = utils.generate_genome_tag(args.reference_path)
    motifs_def = atac_feature_ref.from_motif_list(motif_names, genomes)
    tf_matrix = cr_matrix.CountMatrix(motifs_def, barcodes, tf_peak_matrix * peak_matrix.m)

    # perform MAD-zscoring of proportion values
    propZ_matrix = np.array(tf_matrix.m / peak_matrix.m.sum(axis=0))
    propZ_matrix = MADzscore(propZ_matrix)

    outs.coerce_strings()

    # save to h5 and csv
    tf_matrix.save_h5_file(outs.filtered_tf_bc_matrix, sw_version=martian.get_pipelines_version())
    if not os.path.exists(outs.filtered_tf_bc_matrix_mex):
        os.mkdir(outs.filtered_tf_bc_matrix_mex)
    atac_matrix.save_mex(tf_matrix,
                         outs.filtered_tf_bc_matrix_mex,
                         feature_type=cr_lib_constants.ATACSEQ_LIBRARY_DERIVED_TYPE,
                         sw_version=martian.get_pipelines_version())
    # save propZ matrix as gz
    np.savetxt(outs.tf_propZ_matrix, propZ_matrix)
