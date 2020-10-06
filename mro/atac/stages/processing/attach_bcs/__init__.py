"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Put the sample index barcode and 10X barcode into read tags in a BAM file.
Screen each 10X barcode against the given barcode whitelist
"""

import itertools
import json
import random

import martian
import numpy as np

import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.seq as tk_seq
from barcodes import load_barcode_whitelist, get_read_barcode, stringent_read_filter
from constants import PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG, SELF_FIVE_PRIME_POS_TAG, \
    MATE_FIVE_PRIME_POS_TAG, MATE_MAPPING_QUALITY_TAG, TENX_PRODUCT_NAME
from tenkit.constants import SAMPLE_INDEX_TAG, SAMPLE_INDEX_QUAL_TAG, ILLUMINA_QUAL_OFFSET, TRIM_TAG, TRIM_QUAL_TAG
from tools import open_maybe_gzip, create_bam_infile
from tools.peaks import compute_five_prime_coords
from constants import MAX_I32

__MRO__ = """
stage ATTACH_BCS(
    in  string barcode_whitelist,
    in  bam[]  align,
    in  map[]  chunks,
    in  bool   paired_end,
    in  bool   exclude_non_bc_reads,
    in  float  bc_confidence_threshold,
    in  json   bc_counts,
    out bam[]  output,
    out int    perfect_read_count,
    src py     "stages/processing/attach_bcs",
) split (
    in  bam    align_chunk,
    in  map    chunk,
)
"""

DNA_ALPHABET = 'AGCT'
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}
ALPHABET_MINUS['N'] = set(DNA_ALPHABET)
MAXDIST_CORRECT = 2


class GlobalFivePrimePosTagger:
    """
    Computes the location of the 5' most end of each read in chromosome-stiched global
    reference coordinates including soft-clipped regions. This will be useful in
    marking duplicates accurately.

    Attributes:

    """

    def __init__(self, bam_in):
        """
        Args:
            bam_in(pysam.AlignmentFile): Pysam bam reader object
        """
        chrom_lengths = bam_in.lengths  # Lengths of chromosomes in the same order as pysam.AlignmentFile.references
        self.offsets = [0] * bam_in.nreferences
        for i in xrange(1, bam_in.nreferences):
            self.offsets[i] = self.offsets[i - 1] + chrom_lengths[i - 1]

    def compute_global_coords(self, read):
        """
        For unmapped read the global coordinate is -1
        For mapped reads, compute the 5' position in chromosome coordinates and add the offset

        Args:
            read(pysam.AlignedSegment): Pysam read object
        """
        # BAM tags are unfortunately of type int32, so wrap around
        # Overlaps are okay as samtools sort will use read position as a second index
        # when we ask it to sort by some tag
        if read.is_unmapped:
            return -1
        else:
            return (compute_five_prime_coords(read) + self.offsets[read.reference_id]) % MAX_I32

    def tag_reads(self, reads):
        """
        Adds an integer tag to the read, denoting the global position of the most 5' base
        Args:
            reads(pysam.AlignedSegment): List of Pysam read objects sharing the same qname
        """
        coords = {True: -1, False: -1}
        mapqs = {True: -1, False: -1}
        for read in reads:
            assert read.query_name == reads[0].query_name
            if not read.is_secondary:
                coords[read.is_read1] = self.compute_global_coords(read)
                mapqs[read.is_read1] = read.mapping_quality

        for read in reads:
            read.set_tag(SELF_FIVE_PRIME_POS_TAG, self.compute_global_coords(read))
            read.set_tag(MATE_FIVE_PRIME_POS_TAG, coords[not read.is_read1])
            read.set_tag(MATE_MAPPING_QUALITY_TAG, mapqs[not read.is_read1])


def split(args):
    """ Attach BCS to each chunk of the input files """
    chunk_defs = [{'chunk': fastq_chunk, 'align_chunk': aln_chunk} for
                  (fastq_chunk, aln_chunk) in
                  zip(args.chunks, args.align)]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):
    """ Pass through the BAM files for sorting """
    outs.output = [chunk.output for chunk in chunk_outs]
    outs.perfect_read_count = sum([chunk.perfect_read_count for chunk in chunk_outs])


def main(args, outs):
    """ Attaches barcodes. Attaches raw barcode to RAW_BC tag and filters those to form set of PROCESSES_BARCODES """

    chunk = args.chunk

    bam_in = create_bam_infile(args.align_chunk)

    bam_out, _ = tk_bam.create_bam_outfile(outs.output, None, None, template=bam_in,
                                           pgs=[tk_bam.make_pg_header(martian.get_pipelines_version(),
                                                                      "attach_bcs",
                                                                      TENX_PRODUCT_NAME)])

    gp_tagger = GlobalFivePrimePosTagger(bam_in)

    if args.barcode_whitelist is None or args.bc_counts is None:
        # If there's no whitelist or counts then all high quality BC reads get allowed.
        barcode_whitelist = None
        wl_idxs = None
        bc_dist = None
    else:
        barcode_whitelist = load_barcode_whitelist(args.barcode_whitelist)

        # Load the bc counts for this GEM group
        counts = json.load(open(args.bc_counts, 'r'))
        counts = counts[str(chunk['gem_group'])]['bc_counts']

        # Prior distribution over barcodes, with pseudo-count
        bc_dist = np.array(counts, dtype=np.float) + 1.0
        bc_dist = bc_dist / bc_dist.sum()
        wl_idxs = {bc: idx for (idx, bc) in enumerate(sorted(list(barcode_whitelist)))}

    # set random seed to get deterministic subsampling
    random.seed(0)

    if chunk['barcode'] is not None:
        processed_barcode_iter = get_raw_processed_barcodes(open_maybe_gzip(chunk['barcode']), barcode_whitelist,
                                                            args.bc_confidence_threshold, chunk['gem_group'],
                                                            chunk['barcode_reverse_complement'], wl_idxs, bc_dist)
        require_barcode_for_stringent = True
    else:
        processed_barcode_iter = itertools.repeat(None)
        require_barcode_for_stringent = False

    if chunk['sample_index'] is not None:
        sample_index_iter = tk_fasta.read_generator_fastq(open_maybe_gzip(chunk['sample_index']))
    else:
        sample_index_iter = itertools.repeat(None)

    if chunk['trim'] is not None:
        trim_iter = tk_fasta.read_generator_fastq(open_maybe_gzip(chunk['trim']), paired_end=True)
    else:
        trim_iter = itertools.repeat(None)

    iters = itertools.izip(processed_barcode_iter, sample_index_iter, trim_iter)

    # First read
    try:
        read = bam_in.next()
    except StopIteration:
        read = None

    # Number of perfect reads -- used to compute down-sampling rates in mark_duplicates
    perfect_read_count = 0

    # Due to secondary alignments, we must apply the tags to all
    # reads with the same cluster name.
    for (barcode_info, sample_index_info, trim_info) in iters:
        tags = []
        read_name = None

        if read is None:
            break

        if barcode_info is not None:
            (bc_read_name, raw_bc_seq, processed_bc_seq, raw_bc_qual) = barcode_info
            tags.append((RAW_BARCODE_TAG, raw_bc_seq))
            tags.append((RAW_BARCODE_QUAL_TAG, raw_bc_qual))
            if processed_bc_seq is not None:
                tags.append((PROCESSED_BARCODE_TAG, processed_bc_seq))
            read_name = bc_read_name.split()[0]

        if sample_index_info is not None:
            (si_read_name, seq, qual) = sample_index_info
            tags.append((SAMPLE_INDEX_TAG, seq))
            tags.append((SAMPLE_INDEX_QUAL_TAG, qual))

            if read_name is not None:
                if si_read_name.split()[0] != read_name:
                    martian.log_info("mismatch: si_read_name: %s, bam_read_name: %s" % (si_read_name, read_name))
                assert (si_read_name.split()[0] == read_name)
            else:
                read_name = si_read_name.split()[0]

        r1_tags = tags
        r2_tags = list(r1_tags)

        if trim_info is not None:
            (trim1_read_name, trim1_seq, trim1_qual, trim2_read_name, trim2_seq, trim2_qual) = trim_info
            if len(trim1_seq) > 0:
                r1_tags.append((TRIM_TAG, trim1_seq))
                r1_tags.append((TRIM_QUAL_TAG, trim1_qual))
            if len(trim2_seq) > 0:
                r2_tags.append((TRIM_TAG, trim2_seq))
                r2_tags.append((TRIM_QUAL_TAG, trim2_qual))

        reads_attached = 0
        reads_to_attach = []

        while read.query_name == read_name or read_name is None:
            tags = r1_tags if read.is_read1 else r2_tags
            if len(tags) > 0:
                existing_tags = read.tags
                existing_tags.extend(tags)
                read.tags = existing_tags

            if reads_to_attach and (
                    read.query_name != reads_to_attach[0].query_name or reads_to_attach[0].query_name is None):
                gp_tagger.tag_reads(reads_to_attach)
                reads_attached += len(reads_to_attach)
                for r in reads_to_attach:
                    if stringent_read_filter(r, require_barcode_for_stringent):
                        perfect_read_count += 1

                    if args.exclude_non_bc_reads:
                        if not (get_read_barcode(r) is None):
                            bam_out.write(r)
                    else:
                        bam_out.write(r)
                reads_to_attach = []

            reads_to_attach.append(read)

            try:
                read = bam_in.next()

            except StopIteration:
                read = None
                break

        gp_tagger.tag_reads(reads_to_attach)
        reads_attached += len(reads_to_attach)
        for r in reads_to_attach:
            if stringent_read_filter(r, require_barcode_for_stringent):
                perfect_read_count += 1

            if args.exclude_non_bc_reads:
                if not (get_read_barcode(r) is None):
                    bam_out.write(r)
            else:
                bam_out.write(r)
        # We may have more than 2 reads if there was a
        # secondary alignment, but less than 2 means
        # something went wrong
        assert (reads_attached >= 2)

    outs.perfect_read_count = perfect_read_count
    bam_out.close()


def get_raw_processed_barcodes(barcode_file, barcode_whitelist, bc_confidence_threshold, gem_group,
                               barcodes_reverse_complement, wl_idxs, wl_dist):
    """Stream through the raw barcodes and generate processed barcodes (which may be none)"""
    bc_iterator = tk_fasta.read_generator_fastq(barcode_file)
    for (name, seq, qual) in bc_iterator:
        if barcodes_reverse_complement:
            seq = tk_seq.get_rev_comp(seq)
            qual = qual[::-1]
        if barcode_whitelist is None:
            corrected_bc = None if ('N' in seq) else seq
        else:
            corrected_bc = correct_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist, MAXDIST_CORRECT)
        if corrected_bc is not None:
            corrected_bc = '{}-{}'.format(corrected_bc, gem_group)
        yield (name, seq, corrected_bc, qual)


def gen_nearby_seqs(seq, qvs, wl_idxs, maxdist=3):
    """Generate all sequences with at most maxdist changes from seq that are in a provided whitelist, along with the
    quality values of the bases at the changed positions.
    """
    allowed_indices = [i for i in range(len(seq)) if seq[i] != 'N']
    required_indices = tuple([i for i in range(len(seq)) if seq[i] == 'N'])
    mindist = len(required_indices)
    if mindist > maxdist:
        return

    for dist in range(mindist + 1, maxdist + 1):
        for modified_indices in itertools.combinations(allowed_indices, dist - mindist):
            indices = set(modified_indices + required_indices)
            error_probs = qvs[np.array(list(indices))]
            for substitutions in itertools.product(
                    *[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
                new_seq = ''.join(substitutions)
                if new_seq in wl_idxs:
                    yield new_seq, error_probs.sum()


def correct_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist, maxdist=3):
    """Estimate the correct barcode given an input sequence, base quality scores, a barcode whitelist, and a prior
    distribution of barcodes.  Returns the corrected barcode if the posterior likelihood is above the confidence
    threshold, otherwise None.  Only considers corrected sequences out to a maximum Hamming distance of 3.
    """
    wl_cand = []
    likelihoods = []

    qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET
    # Restrict QVs within a narrow range to prevent one base error from overly dominating others
    qvs[qvs < 3.0] = 3.0
    qvs[qvs > 40.0] = 40.0

    if seq in wl_idxs:
        if (qvs > 24).all():
            return seq

        wl_cand.append(seq)
        likelihoods.append(wl_dist[wl_idxs[seq]])

    for test_str, error_probs in gen_nearby_seqs(seq, qvs, wl_idxs, maxdist):
        idx = wl_idxs.get(test_str)
        p_bc = wl_dist[idx]
        log10p_edit = error_probs / 10.0
        likelihoods.append(p_bc * (10 ** -log10p_edit))
        wl_cand.append(test_str)

    posterior = np.array(likelihoods)
    posterior /= posterior.sum()

    if len(posterior) > 0:
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]
    return None
