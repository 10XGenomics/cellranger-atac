#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import numpy as np

def get_genomes_from_feature_ref(fref):
    return sorted(set(f.tags.get('genome', '') for f in fref.feature_defs))

def get_cell_associated_barcodes(genomes, filtered_barcodes_csv):
    """ Get cell-associated barcodes by genome.
    Args:
      genomes (list of str): Genome names.
      filtered_barcodes_csv (str): Path to CSV file.
    Returns:
      dict of (str, set): Map genome to list of cell-assoc barcodes. Empty-string key is for all genomes."""
    cell_bcs = {}
    for genome in genomes:
        # Get all cell-assoc barcodes (ignoring genome) for the "" (blank) genome string
        cell_bcs[genome] = cr_utils.get_cell_associated_barcode_set(filtered_barcodes_csv,
                                                                    genome)
    # All cell-associated barcodes
    cell_bcs[''] = reduce(lambda x,y: x | y, cell_bcs.itervalues(), set())
    return cell_bcs

def target_counts(max_count, divisions, round_max_count):
    """ Construct a list of target subsampling counts.
    Args:
      max_count (int/float): the largest target count
      divisions (int): desired number of targets, including max_count
      round_max_count (bool) = final target is round(max_count) instead of max_count
    Returns:
      list of mixed int/float: target subsampling counts (possibly fractional)"""
    increment = 1.0 / float(divisions)
    percentiles = np.arange(increment, 1.0, step=increment)
    targets = map(round, max_count * percentiles)
    if round_max_count:
        targets.append(round(max_count))
    else:
        targets.append(max_count)
    return targets


def make_subsamplings(targets, n_cells_per_lib, library_type, lib_indexes,
                      depth_type, usable_count_per_lib, raw_count_per_lib = None):
    """
    Args:
      targets (list of int/float): target subsampling counts
      n_cells_per_lib (list of int): number of cell-associated barcodes per library
      library_type: libraries of this type will be used to make subsamplings
      lib_indexes: indices into *_per_lib corresponding to this library_type
      depth_type: one of cellranger.constants.*_SUBSAMPLE_TYPE
      usable_count_per_lib: usable read count per library
      raw_count_per_lib: raw read count per library (needed for 'raw' depth_type only)
    Returns:
      list of dict: target subsampling parameters for each target count"""
    n_libraries = len(n_cells_per_lib)
    subsamplings = []

    for target_rppc in sorted(set(targets)):
        if depth_type == cr_constants.RAW_SUBSAMPLE_TYPE:
            # Infer the usable depth required to achieve this raw depth
            usable_read_fracs = usable_count_per_lib.astype(float) / raw_count_per_lib
            target_usable_counts = target_rppc * n_cells_per_lib * usable_read_fracs
        else:
            target_usable_counts = target_rppc * n_cells_per_lib

        # Zero out libraries of the other types
        rates = np.zeros(n_libraries, dtype=float)
        rates[lib_indexes] = target_usable_counts[lib_indexes].astype(float) \
            / usable_count_per_lib[lib_indexes]

        # Clamp rates that are close to 1 to 1
        rates[np.absolute(rates - 1) < 1e-3] = 1

        # Zero out the libraries for which we have fewer reads than the target
        rates[rates > 1] = 0.0

        enough_data = np.any((rates > 0) & (rates <= 1))
        if not enough_data:
            rates = np.zeros(len(rates))

        subsamplings.append({
            'library_type': library_type,
            'subsample_type': depth_type,
            'target_read_pairs_per_cell': int(round(target_rppc)),
            'library_subsample_rates': list(map(float, rates)),
        })

    return subsamplings
