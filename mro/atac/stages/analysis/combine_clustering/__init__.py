#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#

import os
import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.io as cr_io

from constants import ALLOWED_FACTORIZATIONS, CLUSTER_FILE_HEAD

__MRO__ = """
stage COMBINE_CLUSTERING(
    in  h5   filtered_matrix,
    in  map  clustering_summary,
    in  path clustered_data,
    in  map  graph_clustering_summary,
    in  path knn_clusters,
    out path clustering,
    out map  clustering_summary,
    src py   "stages/analysis/combine_clustering",
)
"""

def main(args, outs):
    outs.clustering_summary = {}
    if args.filtered_matrix is None:
        outs.clustering = None
        return

    if not os.path.exists(outs.clustering):
        cr_io.mkdir(outs.clustering)

    # NOTE: both graph clustering and normal clustering should have run for given method
    assert args.clustering_summary['h5'].keys() == args.graph_clustering_summary['h5'].keys()

    outs.clustering_summary = {'h5': {}, 'csv': {}}
    for method in args.clustering_summary['h5'].keys():
        if method not in ALLOWED_FACTORIZATIONS:
            raise ValueError("invalid method")
        merge_h5 = [args.clustering_summary['h5'][method], args.graph_clustering_summary['h5'][method]]
        groups = [analysis_constants.ANALYSIS_H5_MAP_CLUSTERING[method], analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP]

        out_method_dir = os.path.join(outs.clustering, method)
        cr_io.mkdir(out_method_dir, allow_existing=True)

        out_clustering_h5 = os.path.join(out_method_dir, "{}_clustering.h5".format(method))
        outs.clustering_summary['h5'][method] = out_clustering_h5
        analysis_io.combine_h5_files(merge_h5, out_clustering_h5, groups)

        _csv1 = os.path.join(args.clustered_data, method, CLUSTER_FILE_HEAD[method] + "_csv")
        _csv2 = os.path.join(args.knn_clusters, method, "clusters_csv")
        out_csv = os.path.join(out_method_dir, method + "_csv")
        cr_io.copytree(_csv1, out_csv, allow_existing=True)
        cr_io.copytree(_csv2, out_csv, allow_existing=True)
        outs.clustering_summary['csv'][method] = out_csv
