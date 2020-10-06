"""
Broadly useful tools for peak detection and transposition site identification.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
import itertools

from constants import SELF_FIVE_PRIME_POS_TAG, MATE_FIVE_PRIME_POS_TAG, MAX_I32, CIGAR_SOFTCLIP
from tenkit.constants import READ_MATE_FAR_DIST


def is_chimeric_fragment(read):
    """Determine whether or not a read comes from a chimeric fragment.  Possible kinds of chimeras:

    - R1/R2 map to different chromosomes
    - R1/R2 map more than 5kb away from each other
    - Both R1 and R2 map to the same strand
    - The earlier read is reverse stranded or the later read is forward stranded
    """
    if read.reference_id != read.next_reference_id:
        return True
    elif abs(read.reference_start - read.next_reference_start) > READ_MATE_FAR_DIST:
        return True
    elif read.is_reverse == read.mate_is_reverse:
        return True
    elif read.reference_start < read.next_reference_start and (read.is_reverse or not read.mate_is_reverse):
        return True
    elif read.next_reference_start < read.reference_start and (read.mate_is_reverse or not read.is_reverse):
        return True
    else:
        return False


def adjusted_position_pairs(read):
    """Given a single read, estimate the position of the Tn5 binding sites for both this read's transposition event
    and its mate's.

    Uses the read tags to get folded positions for the fragment start and end, then adjusts them to estimate the Tn5
    binding site.
    """
    invalid_results = None, None

    if is_chimeric_fragment(read):
        return invalid_results

    self_pos = read.get_tag(SELF_FIVE_PRIME_POS_TAG)
    mate_pos = read.get_tag(MATE_FIVE_PRIME_POS_TAG)

    # Because the tag-based positions are capped, we adjust and re-mod with MAX_I32 to correct fragments
    # that happen to cross over the position roll-over
    if abs(mate_pos - self_pos) > (MAX_I32 // 2):
        self_pos = (self_pos + (MAX_I32 // 2)) % MAX_I32
        mate_pos = (mate_pos + (MAX_I32 // 2)) % MAX_I32

    diff = mate_pos - self_pos

    self_truepos = compute_five_prime_coords(read)
    if read.is_reverse:
        stop = self_truepos
        start = self_truepos + diff
    else:
        start = self_truepos
        stop = self_truepos + diff

    if start < 0:
        return invalid_results

    return start + 4, stop - 5


def compute_five_prime_coords(read):
    """
    Computes the 5' position in chromosome coordinates
    Assumes that the read is mapped and the CIGAR exists.

    Args:
        read(pysam.AlignedSegment): Pysam read object
    """
    cigar = read.cigartuples
    if read.is_reverse:
        # Add the suffix clip to the alignment end position
        suffix_clip = sum([x[1] for x in itertools.takewhile(lambda x: x[0] == CIGAR_SOFTCLIP, reversed(cigar))])
        return read.reference_end + suffix_clip
    else:
        # Subtract the prefix clip from the alignment position
        prefix_clip = sum([x[1] for x in itertools.takewhile(lambda x: x[0] == CIGAR_SOFTCLIP, cigar)])
        return read.reference_start - prefix_clip
