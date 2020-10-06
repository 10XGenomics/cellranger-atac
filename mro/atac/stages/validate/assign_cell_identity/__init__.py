"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Analyze the clusters and assign identity to the clusters
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import cellranger.matrix as cr_matrix
import martian

__MRO__ = """
stage ASSIGN_CELL_IDENTITY(
    in  string   factorization,
    in  int      minclusters,
    in  int      maxclusters,
    in  csv      truth_counts,
    in  h5       peak_matrix,
    in  path     analysis,
    in  h5       annotated_peak_matrix,
    in  txt      cell_annotations
    out map      summary,
    out path     assignments,
    src py       "stages/validate/assign_cell_identity",
)
"""
CLUSTERING_METHODS = {'pca': 'kmeans', 'lsa': 'kmeans', 'plsa': 'kmeans'}

def split(args):
    mem_gb = 12.0
    return {'chunks': [], 'join': {'__mem_gb': mem_gb}}

def pearson(v1, v2):
    return pearsonr(v1, v2)[0]

def spearman(v1, v2):
    return spearmanr(v1, v2)[0]

def bhattacharya_coeff(p, q):
    assert p.shape == q.shape
    p1 = p * 1.0 / p.sum()
    q1 = q * 1.0 / q.sum()
    return np.sum(np.sqrt(p1 * q1))

def correlate_cluster_celltypes(clabels, matrix, truth_counts, truth_cell_types):
    """ Aggregate signal across cells from same cluster and correlate with ground
    truth signal for various cell types. Pick strongest correlation as cell identity
    for the cluster. This works most of the time except when disambiguiating close
    cell types such as CD4 and CD8"""

    labels = np.unique(clabels)
    heatmap = np.zeros((len(labels), len(truth_cell_types)))
    if clabels.size != matrix.shape[0] or matrix.shape[1] != truth_counts.shape[0] \
            or truth_counts.shape[1] != len(truth_cell_types):
        raise ValueError("incompatible dimensions of input data")
    for i, label in enumerate(labels):
        agg = np.squeeze(np.array(matrix[clabels == label, :].todense().sum(axis=0)))
        for j in range(len(truth_cell_types)):
            heatmap[i, j] = spearman(agg, truth_counts[:, j])

    # pick closest type:
    mapping = heatmap.argmax(axis=1)
    nclabels = [truth_cell_types[mapping[label - 1]] for label in clabels]
    return nclabels, heatmap

def reduce_on_columns(matrix, column_labels):
    """Generate a matrix with reduction on columns with identical labels"""

    unique_column_labels = np.unique(column_labels)
    ncells, nsites = matrix.shape
    collapsed_matrix = np.zeros((ncells, unique_column_labels.size))
    for i, label in enumerate(unique_column_labels):
        subselect = (column_labels == label)
        collapsed_matrix[:, i] = np.squeeze(matrix[:, subselect].sum(axis=1))
    return collapsed_matrix * 1.0 / (collapsed_matrix.sum(axis=1, keepdims=True) + 1e-28)

def main(args, outs):
    martian.throw("No chunks defined.")

def join(args, outs, chunk_defs, chunk_outs):
    # IO
    if not os.path.exists(args.analysis):
        raise IOError("analysis files not found")

    # read truth counts file and cell types
    if not os.path.exists(args.truth_counts):
        raise IOError("truth counts not found")
    truth_df = pd.read_csv(args.truth_counts, sep=",")
    truth_counts = truth_df.values
    truth_cell_types = truth_df.columns.values.tolist()
    if len(truth_cell_types) == 0:
        raise ValueError("no ground truth cell types provided")

    # sort peak_matrix on barcodes
    if not os.path.exists(args.peak_matrix):
        raise IOError("peak matrix not found")
    peak_matrix_h5 = cr_matrix.CountMatrix.load_h5_file(args.peak_matrix)
    barcodes = peak_matrix_h5.bcs
    perm_ = range(len(barcodes))
    barcodes, perm_ = zip(*sorted(zip(barcodes, perm_)))
    peak_matrix = peak_matrix_h5.m.T[perm_, :]

    if not os.path.exists(outs.assignments):
        os.mkdir(outs.assignments)

    outs.summary = {}
    outs.summary['celltypes'] = truth_cell_types

    # make epinomics heatmaps
    if os.path.exists(args.annotated_peak_matrix):
        if not os.path.exists(args.cell_annotations):
            raise IOError("cell annotations not provided")
        annotated_peak_matrix_h5 = cr_matrix.CountMatrix.load_h5_file(args.annotated_peak_matrix)
        barcodes2 = annotated_peak_matrix_h5.bcs
        perm_ = range(len(barcodes2))
        _, perm_ = zip(*sorted(zip(barcodes2, perm_)))
        barcodes2 = np.sort(barcodes2)
        annotated_peak_matrix = annotated_peak_matrix_h5.m.T[perm_, :]

        cell_annotations = np.squeeze(pd.read_csv(args.cell_annotations, header=None, sep="\t").values.flatten())
        if cell_annotations.size != annotated_peak_matrix.shape[1]:
            raise ValueError("number of annotations is inconsistent with sites in peak matrix file")

        annotated_heatmap = reduce_on_columns(annotated_peak_matrix, cell_annotations)
        annotated_heatmap_file = os.path.join(outs.assignments, "annotation_map.csv")
        bc = pd.DataFrame({'Barcode': barcodes2})
        annotated_heatmap_df = pd.concat([bc, pd.DataFrame(annotated_heatmap, columns=np.unique(cell_annotations))], axis=1)
        annotated_heatmap_df.to_csv(annotated_heatmap_file, sep=",", index=False)
        outs.summary['annotations'] = annotated_heatmap_file

    method = args.factorization
    if 'assignments' not in outs.summary:
        outs.summary['assignments'] = {}
    tsne_file = os.path.join(args.analysis, "tsne", "2_components", "projection.csv")
    tsne_df = pd.read_csv(tsne_file, sep=",")
    cell_barcodes = tsne_df['Barcode'].values.tolist()
    tsne1 = tsne_df['TSNE-1'].values.tolist()
    tsne2 = tsne_df['TSNE-2'].values.tolist()
    cell_barcodes, tsne1, tsne2 = zip(*sorted(zip(cell_barcodes, tsne1, tsne2)))

    # for each clustering
    # associate cluster to cell_barcodes in tsne projection
    for ncluster in range(args.minclusters, args.maxclusters + 1):
        if ncluster not in outs.summary['assignments']:
            outs.summary['assignments'][ncluster] = {}
        filetag = "{}_{}_{}".format(method, CLUSTERING_METHODS[method], ncluster)
        cluster_file = os.path.join(args.analysis, "clustering",
                                    "{}_{}_clusters".format(CLUSTERING_METHODS[method], ncluster),
                                    "clusters.csv")

        cluster_df = pd.read_csv(cluster_file, sep=",")
        cell_barcodes2 = cluster_df['Barcode'].values.tolist()
        if sorted(cell_barcodes) != sorted(cell_barcodes2):
            raise ValueError("clustering and tsne use different cell_barcodes")
        clabels = cluster_df['Cluster'].values.tolist()

        # canonicalize
        _, clabels = zip(*sorted(zip(cell_barcodes2, clabels)))
        # make heatmaps for pearson corr
        # find index of cell_barcodes in all barcodes
        setbc = set(cell_barcodes)
        subselect = np.array([True if b in setbc else False for b in barcodes])
        nclabels, heatmap = correlate_cluster_celltypes(np.array(clabels), peak_matrix[subselect, :], truth_counts, truth_cell_types)
        unified_file = os.path.join(outs.assignments, "{}_unified.csv".format(filetag))
        unified_df = pd.DataFrame({'Barcode': cell_barcodes, 'x': tsne1, 'y': tsne2,
                                   'Cluster': clabels,
                                   'Celltype': nclabels})
        unified_df.to_csv(unified_file, sep=",", header=True, index=False)

        # write heatmap for each clustering
        heatmap_file = os.path.join(outs.assignments, "{}_heatmap.csv".format(filetag))
        np.savetxt(heatmap_file, heatmap, delimiter=",")
        outs.summary['assignments'][ncluster][method] = [unified_file, heatmap_file]

        assert set(barcodes2) == set(barcodes)

    # graph clustering --
    cluster_key = 'graphclust'
    if cluster_key not in outs.summary['assignments']:
        outs.summary['assignments'][cluster_key] = {}
    filetag = "{}_{}".format(method, 'graphclust')
    cluster_file = os.path.join(args.analysis, "clustering", "graphclust", "clusters.csv")

    cluster_df = pd.read_csv(cluster_file, sep=",")
    cell_barcodes2 = cluster_df['Barcode'].values.tolist()
    if sorted(cell_barcodes) != sorted(cell_barcodes2):
        raise ValueError("clustering and tsne use different cell_barcodes")
    clabels = cluster_df['Cluster'].values.tolist()

    # canonicalize
    _, clabels = zip(*sorted(zip(cell_barcodes2, clabels)))
    # make heatmaps for pearson corr
    # find index of cell_barcodes in all barcodes
    setbc = set(cell_barcodes)
    subselect = np.array([True if b in setbc else False for b in barcodes])
    nclabels, heatmap = correlate_cluster_celltypes(np.array(clabels), peak_matrix[subselect, :], truth_counts, truth_cell_types)

    unified_file = os.path.join(outs.assignments, "{}_unified.csv".format(filetag))
    unified_df = pd.DataFrame({'Barcode': cell_barcodes, 'x': tsne1, 'y': tsne2,
                               'Cluster': clabels,
                               'Celltype': nclabels})
    unified_df.to_csv(unified_file, sep=",", header=True, index=False)

    # write heatmap for each clustering
    heatmap_file = os.path.join(outs.assignments, "{}_heatmap.csv".format(filetag))
    np.savetxt(heatmap_file, heatmap, delimiter=",")
    outs.summary['assignments'][cluster_key][method] = [unified_file, heatmap_file]
