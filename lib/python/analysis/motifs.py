"""
Tools for analyzing motifs in the genome.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

import tempfile
import pyfasta
import MOODS
import MOODS.scan
import MOODS.tools
import MOODS.parsers

from Bio import motifs
from tools import ReferenceManager


ALLOWED_FORMATS = ['bed', 'motif', 'position', 'count', 'binary', 'binary-bed']
class Motifs:

    def __init__(self, ref_path, bg=None):
        ref_manager = ReferenceManager(ref_path)
        self.all_motifs = []
        if ref_manager.motifs is not None:
            with open(ref_manager.motifs, "r") as infile:
                self.all_motifs = list(motifs.parse(infile, "jaspar"))

        # for large sequence header, only keep the text before the first space
        self.genome_seq = pyfasta.Fasta(ref_manager.fasta, key_fn=lambda x: x.split()[0])
        self.bg = bg

    def scan_motif_from_bed(self, peaks_iter, tf_genes=None, out_file=None, out_format='bed', use_genome_bg=True,
                            pseudocount=1.0, pvalue=1e-5, window=7):
        """
        :param peaks_iter: peaks iterator, yielding a tuple (chrom, start, end, strand)
        :param tf_genes: a list of gene symbols for motif search. If None, scan all motifs in the collection
        :param out_file: file name of output. Default is no output file
        :param out_format: output format can be bed, position, count or binary
        :param use_genome_bg: whether to generate background model using the genome seq.
                              If False, sequences from the input bed file will be used for background
        :param pseudocount: MOODS scan pseudocounts for the calculation of log-odds scores from matrices
        :param pvalue: motif hits cutoff
        :param window: the scanning window size of the lookahead filtration algorithm
        :return:
        """

        assert out_format in ALLOWED_FORMATS
        if use_genome_bg:
            self.bg = self.get_reference_bg(self.genome_seq)

        motif = self.get_motif_of(tf_genes)

        # Each TF gets a matrix for the + and for the - strand, and a corresponding threshold
        (matrices, thresholds) = self._prepare_moods_settings(motif, self.bg, pseudocount, pvalue)

        scanner = MOODS.scan.Scanner(window)  # parameter is the window size
        scanner.set_motifs(matrices, self.bg, thresholds)

        if out_file is not None:
            out = open(out_file, 'w')

        maps = []
        for peak_idx, peak in enumerate(peaks_iter):

            bed_coord = {'chr': peak[0],
                         'start': int(peak[1]),
                         'stop': int(peak[2]),
                         'name': peak_idx,
                         'strand': peak[3]}
            seq = self.genome_seq.sequence(bed_coord, one_based=False)

            # seq is of unicode format, need to convert to str
            results = scanner.scan(str(seq))
            parsed = self._parse_scan_results(results, motif, bed_coord, out_format)

            if out_file is not None:
                out.writelines(['\t'.join(map(str, item)) + '\n' for item in parsed])
            else:
                maps.extend(parsed)

        if out_file is not None:
            out.close()
        else:
            return maps

    def get_motif_of(self, tf_genes):
        '''if tf_genes is None, return all motifs in the collection, otherwise find a match'''
        if tf_genes is None:
            return self.all_motifs

        selected = []
        for motif in self.all_motifs:
            if motif.name in tf_genes:
                selected.append(motif)
        return selected

    @staticmethod
    def get_reference_bg(fa):
        b = bytearray()
        for chrom in fa.keys():
            b.extend(bytes(fa[chrom]))

        bg = MOODS.tools.bg_from_sequence_dna(str(b), 1)
        return bg

    def _prepare_moods_settings(self, jaspar_motifs, bg, pseduocount, pvalue=1e-5):
        """Find hits of list of jaspar_motifs in pyfasta object fasta, using the background distribution bg and
        pseudocount, significant to the give pvalue
        """
        matrices = [self._jaspar_to_moods_matrix(j, bg, pseduocount) for j in jaspar_motifs]
        matrices = matrices + [MOODS.tools.reverse_complement(m) for m in matrices]
        thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
        return matrices, thresholds

    @staticmethod
    def _jaspar_to_moods_matrix(jaspar_motif, bg, pseudocount):
        """Convert JASPAR motif into a MOODS log_odds matrix, using the give background distribution & pseudocounts
        - JASPAR 2010 matrix_only format::

          >MA0001.1 AGL3
          A  [ 0  3 79 40 66 48 65 11 65  0 ]
          C  [94 75  4  3  1  2  5  2  3  3 ]
          G  [ 1  0  3  4  1  0  5  3 28 88 ]
          T  [ 2 19 11 50 29 47 22 81  1  6 ]

        - JASPAR 2010-2014 PFMs format::

          >MA0001.1 AGL3
          0       3       79      40      66      48      65      11      65      0
          94      75      4       3       1       2       5       2       3       3
          1       0       3       4       1       0       5       3       28      88
          2       19      11      50      29      47      22      81      1       6

        """
        with tempfile.NamedTemporaryFile() as fn:
            f = open(fn.name, "w")
            for base in 'ACGT':
                line = " ".join(str(x) for x in jaspar_motif.pwm[base])
                f.write(line + "\n")

            f.close()
            return MOODS.parsers.pfm_to_log_odds(fn.name, bg, pseudocount)

    @staticmethod
    def _parse_scan_results(moods_scan_res, motifs, bed_coord, out_format="bed"):
        """ Parse results of MOODS.scan.scan_dna and return/write
            The default input contains one pair of a single motif: forward and reverse strand

            out_format: if "bed" each hit will be the motif position in bed format, with the start and end as peak coordinates
                        if "motif" same as "bed" except the start and end columns are the motif location
                        if "position" each hit will be a list of positions within a peak region (relative to the start position)
                        if "count" each hit will be an integer of the number of occurrences. NO OUTPUT if count == 0
                        if "binary" each hit will be 1 (any occurrence).  NO OUTPUT if count == 0
        """
        assert out_format in ALLOWED_FORMATS

        all_hits = []
        for (motif_idx, hits) in enumerate(moods_scan_res):

            motif = motifs[motif_idx % len(motifs)]
            strand = "-" if motif_idx >= len(motifs) else "+"

            if len(hits) > 0:

                if out_format == "binary":
                    # for binary format we only care about whether len(hits)>0
                    record = [motif_idx % len(motifs), bed_coord['name'], 1 if len(hits) > 0 else 0]
                    all_hits.append(record)
                    continue

                elif out_format == "count":
                    # for count format we just need len(hits)
                    record = [motif_idx % len(motifs), bed_coord['name'], len(hits)]
                    all_hits.append(record)
                    continue

                elif out_format == "binary-bed":
                    record = [bed_coord['chr'], bed_coord['start'], bed_coord['stop'], motif.name]
                    all_hits.append(record)
                    continue

                else:
                    # for bed or position format, we ouput each individual hit
                    # output bed file of [chr start end motifname score strand pos pos+motiflen]

                    for h in hits:
                        if out_format == "bed":
                            motif_start = bed_coord['start'] + int(h.pos)
                            motif_end = bed_coord['start'] + int(h.pos) + motif.length
                            score = round(h.score, 4)
                            record = [bed_coord['chr'], bed_coord['start'], bed_coord['stop'], motif.name, score, strand, motif_start, motif_end]

                        elif out_format == "motif":
                            motif_start = bed_coord['start'] + int(h.pos)
                            motif_end = bed_coord['start'] + int(h.pos) + motif.length
                            score = round(h.score, 4)
                            record = [bed_coord['chr'], motif_start, motif_end, motif.name, score, strand]

                        else:
                            strand = 1 if strand == "-" else -1
                            record = [motif_idx % len(motifs), bed_coord['name'], int(h.pos), round(h.score, 4), strand]

                        all_hits.append(record)

        return all_hits
