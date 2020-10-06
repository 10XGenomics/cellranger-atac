"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Produces a websummary output as a single html page with embedded plots.
"""

from __future__ import division

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams['figure.figsize'] = (9.0, 6.0)
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tools.websummary import HTMLoutput
from metrics import MetricAnnotations

__MRO__ = """
stage CREATE_VALIDATION_PLOTS(
    in  path  assignments,
    in  map   assignment_summary,
    out html  plots,
    src py    "stages/validate/create_validation_plots",
)
"""

def highlight_clusters_on_tsne_plot(metadata, key, assignment_summary, celltypes, annotations):
    """Plot tsne coordinates of barcodes after dimensionality reduction and highlight clusters"""

    # infer from summary
    methods = assignment_summary.keys()
    nmethods = len(methods)
    titles = list(methods)

    fig, axes = plt.subplots(figsize=(6 * nmethods, 36), nrows=6, ncols=nmethods)
    if nmethods == 1:
        axes = np.expand_dims(axes, axis=1)
    cm = plt.get_cmap('hsv')

    annotation_map = pd.read_csv(annotations, sep=",")
    cell_annotations = annotation_map.columns.values.tolist()[1:]
    all_barcodes = annotation_map['Barcode'].values.tolist()
    if all_barcodes != sorted(all_barcodes):
        raise ValueError("wrongly ordered barcodes, should presort them")
    annotations_matrix = annotation_map.drop('Barcode', 1).values

    tsne_results_fixed = None
    method_fixed = None
    lmax_fixed = 0
    lmin_fixed = 0
    # infer from file read in
    col = 0
    for axs, title, method in zip(axes.transpose(), titles, methods):
        unified_file = assignment_summary[method][0]
        heatmap_file = assignment_summary[method][1]

        # results
        unified_df = pd.read_csv(unified_file, sep=",")
        tsne_results = unified_df[['x', 'y']].values
        if col == 0:
            method_fixed = method
            tsne_results_fixed = tsne_results.copy()
            lmin_fixed = min(np.min(tsne_results_fixed[:, 0]), np.min(tsne_results_fixed[:, 1])) - 1
            lmax_fixed = max(np.max(tsne_results_fixed[:, 0]), np.max(tsne_results_fixed[:, 1])) + 1
        col += 1
        barcodes = unified_df['Barcode'].values.tolist()
        clabels = unified_df['Cluster'].values.tolist()
        identity = unified_df['Celltype'].values.tolist()

        # coloring and labeling
        ncluster = len(set(clabels))
        colors = [cm(1. * i / ncluster) for i in range(ncluster)]
        use_colors = [colors[l - 1] for l in clabels]
        ncelltypes = len(celltypes)
        cell_colors = [cm(1. * i / ncelltypes) for i in range(ncelltypes)]
        use_cell_colors = [cell_colors[celltypes.index(c)] for c in identity]

        # plot args
        kwargs = {'lw': 0,
                  'alpha': 0.8,
                  's': 8,
                  'rasterized': True}

        # plot tsne and color clusters
        lmin = min(np.min(tsne_results[:, 0]), np.min(tsne_results[:, 1])) - 1
        lmax = max(np.max(tsne_results[:, 0]), np.max(tsne_results[:, 1])) + 1
        axs[0].scatter(x=tsne_results[:, 0], y=tsne_results[:, 1], color=use_colors, **kwargs)
        # get one representative barcode per group
        for i, group in enumerate(set(use_colors)):
            rep = use_colors.index(group)
            count = use_colors.count(group)
            axs[0].scatter(x=tsne_results[rep, 0], y=tsne_results[rep, 1],
                           c=use_colors[rep], label=str(clabels[rep] - 1) + "(" + str(count) + ")", **kwargs)
        axs[0].legend(loc='lower left', ncol=4, fontsize=8)
        axs[0].set_xlabel('tsne_1')
        axs[0].set_ylabel('tsne_2')
        axs[0].set_title("{}: clusters".format(title))
        axs[0].set_xlim([lmin, lmax])
        axs[0].set_ylim([lmin, lmax])

        # color using identity
        axs[1].scatter(x=tsne_results[:, 0], y=tsne_results[:, 1], color=use_cell_colors, **kwargs)
        # get one representative barcode per group
        for i, group in enumerate(set(use_cell_colors)):
            rep = use_cell_colors.index(group)
            count = use_cell_colors.count(group)
            axs[1].scatter(x=tsne_results[rep, 0], y=tsne_results[rep, 1],
                           c=use_cell_colors[rep], label=identity[rep] + "(" + str(count) + ")", **kwargs)
        axs[1].legend(loc='lower left', ncol=2, fontsize=6)
        axs[1].set_xlabel('tsne_1')
        axs[1].set_ylabel('tsne_2')
        axs[1].set_title("{}: assigned cell-types".format(title))
        axs[1].set_xlim([lmin, lmax])
        axs[1].set_ylim([lmin, lmax])

        # plot heatmap
        heatmap = np.genfromtxt(heatmap_file, delimiter=",")
        hm = axs[2].imshow(heatmap, cmap=plt.cm.Oranges, interpolation='none', aspect='auto')
        axs[2].set_yticklabels(range(ncluster), minor=False)
        axs[2].set_yticks(range(ncluster), minor=False)
        axs[2].set_xticklabels(celltypes, minor=False, fontsize=8, rotation='vertical')
        axs[2].set_xticks(range(ncelltypes), minor=False)
        axs[2].set_title("{}: correlation with purified bulk [Corces et al.]".format(title))
        fig.colorbar(hm, ax=axs[2], cmap=plt.cm.Oranges)

        # plot epinomics map
        # barcodes are already sorted
        setbc = set(barcodes)
        subselect = np.array([True if b in setbc else False for b in all_barcodes])
        annotations_matrix_sub = annotations_matrix[subselect, :]

        # sort cells by cluster labels
        permuted_cell_order = range(len(barcodes))
        ordered_clabels, permuted_cell_order = zip(*sorted(zip(clabels, permuted_cell_order)))
        cluster_grouped_annotations = annotations_matrix_sub[permuted_cell_order, :]
        cell_cuts = [ordered_clabels.index(c) for c in set(clabels)]
        axs[3].imshow(cluster_grouped_annotations, cmap=plt.cm.Greys_r, interpolation='none', aspect='auto')
        for c in sorted(cell_cuts)[1:]:
            axs[3].axhline(y=c, linewidth=1, color='w')
        axs[3].tick_params(axis='y', which='both', left='off', right='off')
        axs[3].set_xticklabels(cell_annotations, minor=False, fontsize=8, rotation='vertical')
        axs[3].set_xticks(range(len(cell_annotations)), minor=False)
        axs[3].set_title("{}: signal at cell-type specific peaks".format(title))

        # plot tsne and color clusters
        axs[4].scatter(x=tsne_results_fixed[:, 0], y=tsne_results_fixed[:, 1], color=use_colors, **kwargs)
        # get one representative barcode per group
        for i, group in enumerate(set(use_colors)):
            rep = use_colors.index(group)
            count = use_colors.count(group)
            axs[4].scatter(x=tsne_results_fixed[rep, 0], y=tsne_results_fixed[rep, 1],
                           c=use_colors[rep], label=str(clabels[rep] - 1) + "(" + str(count) + ")", **kwargs)
        axs[4].legend(loc='lower left', ncol=4, fontsize=8)
        axs[4].set_xlabel('tsne_1')
        axs[4].set_ylabel('tsne_2')
        axs[4].set_title("{}: clusters shown on tsne({})".format(title, method_fixed))
        axs[4].set_xlim([lmin_fixed, lmax_fixed])
        axs[4].set_ylim([lmin_fixed, lmax_fixed])

        # color using identity
        axs[5].scatter(x=tsne_results_fixed[:, 0], y=tsne_results_fixed[:, 1], color=use_cell_colors, **kwargs)
        # get one representative barcode per group
        for i, group in enumerate(set(use_cell_colors)):
            rep = use_cell_colors.index(group)
            count = use_cell_colors.count(group)
            axs[5].scatter(x=tsne_results_fixed[rep, 0], y=tsne_results_fixed[rep, 1],
                           c=use_cell_colors[rep], label=identity[rep] + "(" + str(count) + ")", **kwargs)
        axs[5].legend(loc='lower left', ncol=2, fontsize=6)
        axs[5].set_xlabel('tsne_1')
        axs[5].set_ylabel('tsne_2')
        axs[5].set_title("{}: cell-types shown on tsne({})".format(title, method_fixed))
        axs[5].set_xlim([lmin_fixed, lmax_fixed])
        axs[5].set_ylim([lmin_fixed, lmax_fixed])

    fig.suptitle('Clustering results using multiple methods')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return {'figures': [fig],
            'title': "TSNE plots with clustering (ncluster={c})".format(c=key),
            'text': highlight_clusters_on_tsne_plot.__doc__}

def main(args, outs):
    metadata = MetricAnnotations()

    with HTMLoutput(outs.plots, "") as plotout:
        if args.assignment_summary:
            for key in args.assignment_summary['assignments'].keys():
                plotout.add_results(**highlight_clusters_on_tsne_plot(metadata,
                                    key,
                                    args.assignment_summary['assignments'][key],
                                    args.assignment_summary['celltypes'],
                                    args.assignment_summary['annotations']))
