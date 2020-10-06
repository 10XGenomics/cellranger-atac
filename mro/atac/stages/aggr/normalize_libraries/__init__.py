
"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Perform normalization of libraries
"""

# Imports
import pickle
import json
import os
import pysam
from pybedtools import BedTool

from tools.io import merge_keyed_bed, subsample_fragments
from tools import ReferenceManager
import cellranger.io as cr_io
import martian
from constants import PEAK_MERGE_DISTANCE
from utils import get_cell_barcodes

# MRO docstring
__MRO__ = """
stage NORMALIZE_LIBRARIES(
    in  pickle     library_info,
    in  string     reference_path,
    out tsv.gz     fragments,
    out tsv.gz.tbi fragments_index,
    out csv        cell_barcodes,
    out bool       skip_peakcalling,
    out bed        peaks,
    out json       normalization_metrics,
    src py         "stages/aggr/normalize_libraries",
) split (
    in  int        n,
)
"""

def split(args):
    """Each chunk gets a library to downsample, along with normalization parameters"""
    with open(args.library_info, 'r') as f:
        library_info = pickle.load(f)
    return {'chunks': [{'n': group, '__mem_gb': 12} for group in library_info], 'join': {'__mem_gb': 12}}

def main(args, outs):
    """Downsample each fragments file to produce a sorted file, while computing the pre and post complexity metrics"""
    with open(args.library_info, 'r') as f:
        library_info = pickle.load(f)[args.n]

    # read cells
    cell_barcodes = get_cell_barcodes(library_info['cells'], args.reference_path)

    # get chrom key from fasta index
    chrom_order = {}
    ctg_mgr = ReferenceManager(args.reference_path)
    with open(ctg_mgr.fasta_index, 'r') as f:
        for en, line in enumerate(f):
            chrom = line.split('\t')[0]
            chrom_order[chrom] = en

    downsampling_metrics = subsample_fragments(infile=library_info['fragments'],
                                               rate=library_info['rate'],
                                               outfile=outs.fragments,
                                               group=args.n,
                                               cells=cell_barcodes,
                                               kind=library_info['kind'],
                                               key=chrom_order)

    with open(outs.normalization_metrics, 'w') as f:
        json.dump(downsampling_metrics, f, indent=4)

def join(args, outs, chunk_defs, chunk_outs):
    """Merge sorted, downsampled fragments from each chunk,
    emit pre and post normalization sensitivity metrics per library
    and merge input peaks if provided"""

    with open(args.library_info, 'r') as f:
        library_info = pickle.load(f)

    ctg_mgr = ReferenceManager(args.reference_path)

    # Merge cell_barcodes
    cell_barcodes = {}
    for group in library_info:
        cell_barcodes_group = get_cell_barcodes(library_info[group]['cells'], args.reference_path, with_species=True)
        group_suffix = "-{}".format(group)
        for species in cell_barcodes_group.keys():
            if species not in cell_barcodes.keys():
                cell_barcodes[species] = set()
            cell_barcodes[species].update({bc.split("-")[0] + group_suffix for bc in cell_barcodes_group[species]})
    with open(outs.cell_barcodes, 'w') as f:
        for species in cell_barcodes:
            f.write(species + "," + ",".join(cell_barcodes[species]) + "\n")

    # Merge peaks if provided
    input_peaks = [library_info[group]['peaks'] for group in library_info if 'peaks' in library_info[group]]
    if len(input_peaks) == 1:
        cr_io.copy(input_peaks[0], outs.peaks)
        outs.skip_peakcalling = True
    if len(input_peaks) == 0:
        outs.peaks = 0
        outs.skip_peakcalling = False
    if len(input_peaks) > 1:
        outs.skip_peakcalling = True
        # cat
        with open(outs.peaks, 'w') as outf:
            for ip in input_peaks:
                with open(ip, 'r') as inf:
                    for line in inf:
                        outf.write(line)
        # sort
        peaks = BedTool(outs.peaks)
        peaks = peaks.sort(faidx=ctg_mgr.fasta_index)

        # merge
        peaks = peaks.merge(d=PEAK_MERGE_DISTANCE)
        peaks.saveas(outs.peaks)

    # override library name when aggring 1 library:
    if len(library_info) == 1:
        library_info[1]['library_info'] = ""

    # merge the metrics
    normalization_metrics = {}
    for cdef, cout in zip(chunk_defs, chunk_outs):
        with open(cout.normalization_metrics, 'r') as f:
            chunk_metrics = json.load(f)
            for key in chunk_metrics:
                normalization_metrics["{}_Library_{}".format(key, library_info[cdef.n]['library_id'])] = chunk_metrics[key]
                # aggregate some metrics across all libraries
                if key in ['total_pre_normalization', 'total_post_normalization']:
                    if key not in normalization_metrics:
                        normalization_metrics[key] = 0
                    normalization_metrics[key] += chunk_metrics[key]
    with open(outs.normalization_metrics, 'w') as f:
        json.dump(normalization_metrics, f, indent=4)

    # merge the fragments
    base_file, extension = os.path.splitext(outs.fragments)
    if not extension == '.gz':
        raise ValueError('Expecting compressed file output')
    input_tsvs = [str(chunk.fragments) for chunk in chunk_outs]
    merge_keyed_bed(input_tsvs, base_file, threads=martian.get_threads_allocation())

    # index the fragments
    if os.path.getsize(base_file) == 0:
        outs.fragments = None
        outs.fragments_index = None
    else:
        # N.B. tabix_index will automatically compress the input file, adding the .gz suffix
        pysam.tabix_index(base_file, preset='bed', index=outs.fragments_index)
