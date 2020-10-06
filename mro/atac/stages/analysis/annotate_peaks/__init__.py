u'''
 Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

 Annotate peaks with near by genes, TSS/distal position and motif positions.
'''

from __future__ import division

import martian
import utils
from pybedtools import BedTool, create_interval_from_list
from tools import ReferenceManager
from tools.io import combine_csv
import pandas as pd
from constants import TRANSCRIPT_ANNOTATION_GENE_TYPES

__MRO__ = """
stage ANNOTATE_PEAKS(
    in  bed     peaks,
    in  string  reference_path,
    out tsv     peak_annotation,
) split (
    in  int    chunk_start,
    in  int    chunk_end,
) using (
    mem_gb = 5,
)
"""

DISTANCE_LEFT_OF_TSS = 1000
DISTANCE_RIGHT_OF_TSS = 100
DISTANCE_TO_INTERGENIC = 200000

def get_peak_nearby_genes(peak):
    """
    parse peaks with closest TSS and overlaped transcript
    :param peak: can be 6 columns when a peak has closest TSS and also overlaps with a transcript:
                 chr1	937133	937621	ISG15	-616	ISG15
                 or just 5 column if the peak doesn't overlap with any transcript
                 chr1	937133	937621	ISG15	-616
                 the 4th and 6th (if present) column is a string in which gene symbols are separated by comma
                 the 5th column is distance values separated by comma
                 column 1 to 5 are required and column 6 is optional
    :return: Interval object created from a list like ['chr1', '6051602', '6053638', 'KCNAB2;NPHP4', '0,0', 'promoter;promoter']
    """
    distances = [c for c in peak[4].split(',')]
    genes = [c for c in peak[3].split(',')]
    assert len(distances) == len(genes)
    peak_types = []

    # call promoter peaks first
    is_promoter = dict()
    for i, d in enumerate(distances):
        d = int(d)
        if -d <= DISTANCE_LEFT_OF_TSS and d <= DISTANCE_RIGHT_OF_TSS:
            is_promoter[genes[i]] = True

    # call distal peaks
    is_distal = dict()
    for i, di in enumerate(distances):
        d = int(di)
        if -d <= DISTANCE_LEFT_OF_TSS and d <= DISTANCE_RIGHT_OF_TSS:
            peak_types.append('promoter')

        elif abs(d) <= DISTANCE_TO_INTERGENIC:
            # distal peaks, if this peak is already a promoter of the gene being tested, do not annotate it again as a distal peak
            if genes[i] in is_promoter.keys():
                genes[i] = ''
                distances[i] = ''
                peak_types.append('')
            else:
                peak_types.append('distal')
                is_distal[genes[i]] = True
        else:
            genes[i] = ''
            distances[i] = ''
            peak_types.append('')

        if genes.count('') == len(genes):
            genes, distances, peak_types = [], [], []

    # if a peak has an overlapping gene AND it has not been annotated as the promoter peak of that gene , call distal peaks
    if len(list(peak)) > 5:
        for c in peak[5].split(','):
            if c not in is_promoter.keys() and c not in is_distal.keys():
                distances.append('0')
                peak_types.append('distal')
                genes.append(c)

    # if the peak still has not been annotated, it is an intergenic peak
    if peak_types == []:
        peak_types.append('intergenic')
        genes = ['']
        distances = ['']

    # unify
    assert len(genes) == len(peak_types)
    genes, distances, peak_types = zip(*sorted(set(zip(genes, distances, peak_types))))

    return create_interval_from_list([peak[0], peak[1], peak[2], ';'.join(genes), ';'.join(distances), ';'.join(peak_types)])

def annotate_peaks(peaks, ref_path):
    """
    peak to gene annotation strategy:
        1. if a peak overlaps with promoter region (-1kb, + 100) of any TSS, call it a promoter peak
        2. if a peak is within 200kb of the closest TSS, AND if it is not a promoter peak, call it a distal peak
        3. if a peak overlaps of a transcript, AND it is not a promoter nor a distal peak of the gene, call it a distal peak
            This step is optional
        4. call it an intergenic peak
    """

    ref_mgr = ReferenceManager(ref_path)
    tss = BedTool(ref_mgr.tss_track)

    # if tss.bed contains the 7th column (gene type), then apply filter. Otherwise use all tss sites
    if tss.field_count() == 7:
        tss_filtered = tss.filter(lambda x: x[6] in TRANSCRIPT_ANNOTATION_GENE_TYPES).saveas()
    else:
        df_tss = tss.to_dataframe()
        df_tss['gene_type'] = '.'
        tss_filtered = BedTool.from_dataframe(df_tss).saveas()

    # including transcripts.bed is optional
    if ref_mgr.transcripts_track is None:
        transcripts_filtered = BedTool([])
    else:
        transcripts = BedTool(ref_mgr.transcripts_track)
        if transcripts.field_count() == 7:
            transcripts_filtered = transcripts.filter(lambda x: x[6] in TRANSCRIPT_ANNOTATION_GENE_TYPES).saveas()
        else:
            df_tx = transcripts.to_dataframe()
            df_tx['gene_type'] = '.'
            transcripts_filtered = BedTool.from_dataframe(df_tx).saveas()

    # run bedtools closest for peaks against filtered tss, group by peaks and summarize annotations from select columns
    peaks_nearby_tss = peaks.closest(tss_filtered, D='b', g=ref_mgr.fasta_index).groupby(g=[1, 2, 3], c=[7, 11], o=['collapse']).saveas()

    results = []
    peaks_nearby_tss_butno_tx = peaks_nearby_tss.intersect(transcripts_filtered, v=True).saveas()

    # avoid error when no peaks overlap with any transcipts
    if len(peaks_nearby_tss_butno_tx) < len(peaks_nearby_tss):
        peaks_nearby_tss_and_tx = peaks_nearby_tss \
            .intersect(transcripts_filtered, wa=True, wb=True) \
            .groupby(g=[1, 2, 3, 4, 5], c=[9], o=['distinct'])

        for peak in peaks_nearby_tss_and_tx:
            results.append(get_peak_nearby_genes(peak))

    for peak in peaks_nearby_tss_butno_tx:
        results.append(get_peak_nearby_genes(peak))

    return results


def split(args):
    """Compute base background in split and use it in each chunk
    """

    n_peaks = utils.quick_line_count(args.peaks) if args.peaks else 0
    ref_mgr = ReferenceManager(args.reference_path)
    if len(ref_mgr.list_species()) > 1 or n_peaks == 0 or ref_mgr.tss_track is None:
        chunk_def = [{'skip': True}]
        return {'chunks': chunk_def}

    # write rows of each chunk to a new peak file
    mem_in_gb = 4.0
    chunk_def = [{'__mem_gb': mem_in_gb,
                  'skip': False,
                  'chunk_start': chunk[0],
                  'chunk_end': chunk[1]} for chunk in utils.get_chunks(n_peaks, chunks=20)]
    return {'chunks': chunk_def}


def main(args, outs):
    if args.skip or args.chunk_start == args.chunk_end:
        outs.peak_annotation = None
        return

    # compute nearby genes
    peaks = BedTool(args.peaks)
    peak_chunk = BedTool(peaks[args.chunk_start:args.chunk_end]).saveas()
    annotations = annotate_peaks(peak_chunk, args.reference_path)
    if annotations == []:
        outs.peak_annotation = None
    else:
        ref_mgr = ReferenceManager(args.reference_path)
        annotation_bed = BedTool(annotations).sort(g=ref_mgr.fasta_index).saveas()
        with open(outs.peak_annotation, 'w') as out:
            out.write('\t'.join(['peak', 'gene', 'distance', 'peak_type']) + '\n')
            for row in annotation_bed:
                peak = '_'.join(row[:3])
                out.write('\t'.join([peak, row[3], row[4], row[5]]) + '\n')


def join(args, outs, chunk_defs, chunk_outs):
    if not chunk_defs or chunk_defs[0].skip:
        martian.log_info('Skipping peak annotation')
        outs.peak_annotation = None
        return

    chunk_peak_annotations = [chunk.peak_annotation for chunk in chunk_outs if chunk.peak_annotation is not None]
    combine_csv(chunk_peak_annotations, outs.peak_annotation)
