#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

""" ATAC-specific matrix functionality """

import cellranger.atac.feature_ref as atac_feature_ref
import cellranger.library_constants as lib_constants

def save_mex(matrix, base_dir, feature_type, sw_version, compress=False):
    """ Save an ATAC matrix in Matrix Market Exchange format
    Args:
      matrix (CountMatrix): Matrix to write.
      base_dir (str): Path to output directory.
      sw_version (str): Version of this software.
      feature_type: selects the method of saving features based on feature type
    """
    mex_metadata = {
        'software_version': sw_version,
    }

    if feature_type == lib_constants.ATACSEQ_LIBRARY_TYPE:
        save_features = atac_feature_ref.save_features_bed
    elif feature_type == lib_constants.ATACSEQ_LIBRARY_DERIVED_TYPE:
        save_features = atac_feature_ref.save_motifs_tsv
    else:
        raise ValueError('unidentifiable library type')

    matrix.save_mex(base_dir,
                    save_features,
                    metadata=mex_metadata,
                    compress=compress)
