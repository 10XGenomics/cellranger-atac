"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division
import os

import h5py as h5
import martian
import tables
from scipy.sparse import vstack
import pandas as pd
import numpy as np

import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io
import cellranger.matrix as cr_matrix
from cellranger.feature_ref import FeatureReference
from constants import DEFAULT_FACTORIZATION

__MRO__ = """
stage SUMMARIZE_ANALYSIS(
    in  tsv   peak_annotation,
    in  h5    filtered_peak_bc_matrix,
    in  h5    filtered_tf_bc_matrix,
    in  gz    tf_propZ_matrix,
    in  path  reduced_data,
    in  path  reduction_summary,
    in  path  clustering,
    in  map   clustering_summary,
    in  path  tsne_projection,
    in  map   tsne_summary,
    in  path  enrichment_analysis,
    in  map   enrichment_analysis_summary,
    out h5    analysis,
    out path  analysis_csv,
    out h5    feature_bc_matrix,
    src py    "stages/analysis/summarize_analysis",
) split (
)
"""

def split(args):
    # if no cells called
    if args.filtered_peak_bc_matrix is None:
        return {'chunks': [], 'join': {}}

    # peak-bc matrix
    peak_matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_peak_bc_matrix)

    # tf-bc matrix, propZ matrix both are dense
    # As former is stored sparse and latter as dense, the multiplying factor is (3 + 1) = 4
    tf_matrices_mem_gb = 0
    if args.filtered_tf_bc_matrix is not None:
        ntfs, nbcs, nnz = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_tf_bc_matrix)
        BYTES_PER_FLOAT = np.dtype(float).itemsize
        BYTES_PER_GB = 1024**3
        tf_matrices_mem_gb = 4 * BYTES_PER_FLOAT * ntfs * nbcs / BYTES_PER_GB

    return {'chunks': [], 'join': {'__mem_gb': max(int(np.ceil(peak_matrix_mem_gb + tf_matrices_mem_gb)), h5_constants.MIN_MEM_GB)}}

def main(args, outs):
    martian.throw('No chunks defined.')

def join(args, outs, chunk_defs, chunk_outs):
    if args.filtered_peak_bc_matrix is None or not args.reduction_summary['h5'].keys():
        outs.analysis = None
        outs.analysis_csv = None
        outs.feature_bc_matrix = None
        return

    # Make the FBM
    # build joint Peak + TF count matrix for single genomes
    # combine peak annotations for single genome analysis
    peak_annotation = None
    if args.peak_annotation:
        annotations = pd.read_csv(args.peak_annotation, sep='\t')[['gene', 'peak_type']]
        annotations = annotations.replace(np.nan, '', regex=True)
        annotations = annotations.values.astype(str).tolist()
        peak_annotation = []
        for row in annotations:
            genes = row[0].split(";")
            annotation = row[1].split(";")
            promoter = []
            nearby_gene = []
            assert len(annotation) == len(genes)
            for en, kind in enumerate(annotation):
                if kind == 'promoter':
                    promoter += [genes[en]]
                nearby_gene += [genes[en]]
            peak_annotation += [[';'.join(promoter), ';'.join(nearby_gene)]]
    fbm = cr_matrix.CountMatrix.load_h5_file(args.filtered_peak_bc_matrix)
    mapping = None
    if args.filtered_tf_bc_matrix:
        # combine matrices, ensure the barcodes are same and ordered the same way
        tf_matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_tf_bc_matrix)
        assert (fbm.bcs == tf_matrix.bcs).all()
        if peak_annotation is not None:
            fbm.feature_ref = FeatureReference.addtags(fbm.feature_ref, ['promoter', 'nearby_gene'], peak_annotation)
            tf_matrix.feature_ref = FeatureReference.addtags(tf_matrix.feature_ref, ['promoter', 'nearby_gene'])
        combined_feature_defs = FeatureReference.join(fbm.feature_ref, tf_matrix.feature_ref)
        combined_matrix = vstack([fbm.m, tf_matrix.m])
        # explicit map linking rows in diffexp to combined matrix
        mapping = np.zeros((tf_matrix.features_dim, 2))
        for x in range(tf_matrix.features_dim):
            mapping[x, 0] = x
            mapping[x, 1] = x + fbm.features_dim
        fbm = cr_matrix.CountMatrix(combined_feature_defs, fbm.bcs, combined_matrix)
    fbm.save_h5_file(outs.feature_bc_matrix, sw_version=martian.get_pipelines_version())

    # Pytables doesn't support variable len strings, so use h5py first
    with h5.File(outs.feature_bc_matrix, 'r') as matrix, \
            h5.File(outs.analysis, 'w') as out:
        # TODO: copy the first group; fixme when we have a key
        name = matrix.keys()[0]
        matrix.copy(matrix[name], out, name='matrix')

    factorizations = args.reduction_summary['h5'].keys()
    USE_FACTORIZATION = DEFAULT_FACTORIZATION if DEFAULT_FACTORIZATION in factorizations else factorizations[0]
    with tables.open_file(outs.analysis, 'a') as out:
        for summary, key in zip([args.reduction_summary, args.clustering_summary, args.tsne_summary,
                                 args.enrichment_analysis_summary],
                                [USE_FACTORIZATION, 'clustering', 'tsne', 'enrichment']):
            if summary is None or not summary:
                continue
            print(key, summary)
            data_h5 = summary['h5'][USE_FACTORIZATION]
            with tables.open_file(data_h5, 'r') as indata:
                indata.copy_children(indata.root, out.root, recursive=True)
            dirname = os.path.join(outs.analysis_csv, key)
            cr_io.copytree(summary['csv'][USE_FACTORIZATION], dirname)

    # if mapping is present (single genome case), so is the coloring matrix
    if mapping is not None:
        with h5.File(outs.analysis, 'a') as out:
            out.create_dataset('feature_DE_map', data=mapping)
        args.coerce_strings()
        tf_propZ_matrix = np.loadtxt(args.tf_propZ_matrix)
        with h5.File(outs.analysis, 'a') as out:
            out.create_dataset('diffexp_coloring_matrix', data=tf_propZ_matrix)
