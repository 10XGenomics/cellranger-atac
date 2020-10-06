"""
Barcode handling functions, including query functions that override tenkit query functions
(Due to changes in tags used to store barcodes)

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

import os

from constants import RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG, PROCESSED_BARCODE_TAG
from tools.io import open_maybe_gzip
from tenkit.constants import HIGH_CONF_MAPQ, MIN_MATE_OFFSET_DUP_FILTER

def stringent_read_filter(bam_read, require_barcode):
    """Test for a very high-quality read. Only reads satisfying this predicate are
    used when computing summary dup rates to avoid spurious duplicates.
    Reads must have a high MAPQ, a perfect cigar, with a fragment size somewhat
    longer than the read length """

    high_mapq = bam_read.mapq >= HIGH_CONF_MAPQ
    is_primary = not bam_read.is_secondary
    valid_barcode = (not require_barcode) or get_read_barcode(bam_read) is not None
    perfect_cigar = len(bam_read.cigar) == 1
    mate_diff_contig = bam_read.tid != bam_read.rnext
    long_fragment = mate_diff_contig or abs(bam_read.pos - bam_read.mpos) >= MIN_MATE_OFFSET_DUP_FILTER

    return high_mapq and is_primary and valid_barcode and perfect_cigar and long_fragment

def whitelist_mem_gb(fn):
    """Return memory consumption in loading the whitelist file"""
    path = barcode_whitelist_path(fn)
    return os.path.getsize(path) / 1e9

def barcode_whitelist_path(fn):
    """Barcode whitelist is just a text file of valid barcodes, one per line.
    Lines containing the '#' character are ignored"""
    if fn is None:
        return None
    code_path = os.path.dirname(os.path.abspath(__file__))
    if os.path.exists(fn):
        return fn
    elif os.path.exists(os.path.join(code_path, "{}.txt".format(fn))):
        return os.path.join(code_path, "{}.txt".format(fn))
    elif os.path.exists(os.path.join(code_path, "{}.txt.gz".format(fn))):
        return os.path.join(code_path, "{}.txt.gz".format(fn))
    else:
        raise ValueError("Whitelist {} not found".format(fn))


def load_barcode_whitelist(fn, ordered=False):
    """ Barcode whitelists are text files of valid barcodes, one per line.
    Lines containing the '#' character are ignored
    """
    full_path = barcode_whitelist_path(fn)
    if full_path is None:
        return None
    with open_maybe_gzip(full_path, 'r') as infile:
        if ordered:
            barcodes = [line.strip() for line in infile if '#' not in line]
        else:
            barcodes = {line.strip() for line in infile if '#' not in line}
    return barcodes


def get_barcode_gem_group(barcode):
    split = barcode.split("-")
    if len(split) == 2:
        barcode, gem_group = split
    else:
        gem_group = None
    return gem_group


def split_barcode(barcode, return_gg=False):
    """Splits a barcode sequence into part A, part C, and part B sequences."""
    split = barcode.split("-")
    if len(split) == 2:
        # Remove the gem group from the barcode
        barcode, gem_group = split
    else:
        gem_group = None
    assert len(barcode) == 16
    part_a, part_c, part_b = barcode[:7], barcode[7:9], barcode[9:9 + 7]
    if return_gg:
        return part_a, part_c, part_b, gem_group
    return part_a, part_c, part_b


def merge_barcode(part_a, part_c, part_b, gem_group=None):
    barcode = "{}{}{}".format(part_a, part_c, part_b)
    if gem_group is not None:
        return "{}-{}".format(barcode, gem_group)
    return barcode

def query_barcode_subsequences(valid_barcodes):
    """Breaks down a list of valid barcode sequences into unique part A, B, and C
    subsequences, as well as possible gem groups."""
    part_a_seqs = {}
    part_b_seqs = {}
    part_c_seqs = set()
    gem_group_seqs = set()

    for barcode in valid_barcodes:
        part_a, part_c, part_b, gem_group = split_barcode(barcode, return_gg=True)
        part_c_seqs.add(part_c)
        gem_group_seqs.add(gem_group)
        if part_c not in part_a_seqs:
            part_a_seqs[part_c] = set()
        if part_c not in part_b_seqs:
            part_b_seqs[part_c] = set()
        part_a_seqs[part_c].add(part_a)
        part_b_seqs[part_c].add(part_b)

    part_c_seqs = sorted(part_c_seqs)
    for part_c in part_c_seqs:
        part_a_seqs[part_c] = sorted(part_a_seqs[part_c])
        part_b_seqs[part_c] = sorted(part_b_seqs[part_c])

    return part_a_seqs, part_c_seqs, part_b_seqs, gem_group_seqs

def split_barcode_and_gem_group(barcode):
    """Splits a barcode sequence into a barcode and a gem group."""
    split = barcode.split("-")
    if len(split) == 2:
        # Remove the gem group from the barcode
        barcode, gem_group = split
    else:
        gem_group = None
    assert len(barcode) == 16
    return barcode, gem_group

def merge_barcode_and_gem_group(barcode, gem_group):
    """Merges sequence with gem group"""
    if gem_group is not None:
        return "{}-{}".format(barcode, gem_group)
    return barcode

def query_barcodes_and_gem_groups(valid_barcodes):
    """Breaks down a list of valid barcode sequences into barcodes - gem group
    subsequences."""
    barcode_seqs = {}
    gem_group_seqs = set()

    for barcode in valid_barcodes:
        assert '-' in barcode
        barcode_seq, gem_group = barcode.split('-')
        gem_group_seqs.add(gem_group)
        if gem_group not in barcode_seqs:
            barcode_seqs[gem_group] = set()
        barcode_seqs[gem_group].add(barcode_seq)

    gem_group_seqs = sorted(gem_group_seqs)
    for gem_group in gem_group_seqs:
        barcode_seqs[gem_group] = sorted(barcode_seqs[gem_group])

    return barcode_seqs, gem_group_seqs

def _get_tag_if_exists_else_default(read, tag, default=None):
    return read.get_tag(tag) if read.has_tag(tag) else default


def get_read_raw_barcode(read, default=None):
    return _get_tag_if_exists_else_default(read, RAW_BARCODE_TAG, default)


def get_read_barcode_qual(read, default=None):
    return _get_tag_if_exists_else_default(read, RAW_BARCODE_QUAL_TAG, default)


def get_read_barcode(read, default=None):
    tag = _get_tag_if_exists_else_default(read, PROCESSED_BARCODE_TAG, default)
    return default if tag == '' else tag
