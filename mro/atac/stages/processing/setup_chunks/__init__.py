#!/usr/bin/env python
"""Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Chunks up the input fastq data into sets of matched R1, R2, SI, and BC fastq files.

Input args:
sample_id - Arbitrary string labeling the sample
input_mode - one of BCL_PROCESSOR or ILMN_BCL2FASTQ.  BCL_PROCESSOR uses interleaved reads (RA) and labels fastqs
             by sample index.  ILMN_BCL2FASTQ uses R1/R2 and labels by sample name.
downsample - Can be null (for no downsampling) or a map that specifies a subsampling rate < 1.0 (as a float),
             or gigabases of output data requested (as an integer)
sample_def - A list of maps, each with information for each sample to be processed, including:
                bc_in_read - read in which to look for the barcode
                bc_length - length of barcode in read
                read_path - path to FASTQ files
                lanes -
                gem_group -
                unbarcoded - defaults to False if not present, drops barcode information
                library_id
                sample_indices - only for BCL_PROCESSOR, list of valid sample indices, if contains "all", take all
                sample_names - only for ILMN_BCL2FASTQ, list of sample_names to pass on
                subsample_rate
"""

from __future__ import division

import gzip
import itertools
import os.path
import re
import struct

import martian
import tenkit.fasta as tk_fasta
import tenkit.preflight as tk_preflight
import tenkit.safe_json
from constants import BCL2FASTQ_SEQNAMES

__MRO__ = """
stage SETUP_CHUNKS(
    in  string    sample_id          "id of the sample",
    in  map[]     sample_def         "list of dictionary specifying input data",
    in  string    input_mode         "configuration of the input fastqs",
    in  map       downsample         "map specifies either subsample_rate (float) or gigabases (int)",
    out map[]     chunks             "map has barcode, barcode_reverse_complement, sample_index, read1, read2, gem_group, and read_group fields",
    out string[]  read_groups        "list of strings representing read groups",
    out json      downsample_info    "info about downsampling",
    src py        "stages/processing/setup_chunks",
)
"""

def main(args, outs):
    """Combine reads from multiple input FASTQ files, and potentially trim.
       Demultiplex outputs a series of FASTQ files with filenames of the form:
       read-[RA|I1|I2]_si-AGTAACGT_lane-001_chunk_001.fastq[.gz].
    """

    validate_input(args)

    global_subsample_rate = args.downsample.get('subsample_rate', 1.0) if args.downsample is not None else 1.0

    # Predicted input bases
    total_seq_bases = 0

    chunks = []
    read_groups = set()

    for read_chunk in args.sample_def:
        subsample_rate = global_subsample_rate * read_chunk.get('subsample_rate', 1.0)

        bc_in_read = {}
        if read_chunk.get('bc_in_read', None) is not None:
            bc_in_read['bc_in_read'] = read_chunk['bc_in_read']
            bc_in_read['bc_length'] = read_chunk['bc_length']

        path = read_chunk['read_path']
        lanes = read_chunk['lanes']
        gem_group = read_chunk['gem_group']
        unbarcoded = read_chunk.get('unbarcoded', False)
        if unbarcoded:
            martian.log_info('Flagged as unbarcoded: processing as bulk data')

        sample_id = args.sample_id
        library_id = read_chunk.get('library_id', 'MissingLibrary')

        # split on BCL_PROCESSOR / ILMN_BCL2FASTQ
        # the main difference is that BCL_PROCESSOR uses interleaved reads and labels FASTQs by sample index;
        # whereas ILMN_BCL2FASTQ uses R1/R2 and labels by sample name

        if args.input_mode == "BCL_PROCESSOR":
            sample_index_strings, msg = tk_preflight.check_sample_indices(read_chunk)
            if sample_index_strings is None:
                martian.exit(msg)

            sample_seq_bases = 0
            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            for sample_index in sample_index_strings:
                read_paths = find_func(path, "RA", sample_index, lanes)
                for read in read_paths:
                    _, predicted_seq_bases = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases

            martian.log_info("Input data: Predict %f GB from %s" % (sample_seq_bases / 1e9, path))
            total_seq_bases += sample_seq_bases

            for sample_index in sample_index_strings:
                read_paths = find_func(path, "RA", sample_index, lanes)
                # cell barcodes and sample indices are embedded in the index reads
                si_read, bc_read = ("I1", "I2")

                # allow empty sample index case if all reads in lane are same sample
                sis = find_func(path, si_read, sample_index, lanes)
                if sis is None or len(sis) == 0:
                    sis = [None] * len(read_paths)

                barcodes = find_func(path, bc_read, sample_index, lanes)
                if unbarcoded or len(barcodes) == 0:
                    barcodes = [None] * len(read_paths)

                # calculate chunks
                for r, b, si in zip(read_paths, barcodes, sis):
                    (flowcell, lane) = get_run_data(r)
                    if sample_id is not None:
                        rg_string = ':'.join(str(item) for item in [sample_id, library_id, gem_group, flowcell, lane])
                    else:
                        rg_string = 'None:None:None:None:None'
                    new_chunk = {
                        'read1': r, 'read2': None, 'reads_interleaved': True, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

        elif args.input_mode == "ILMN_BCL2FASTQ":
            r1_read, r2_read, si_read, bc_read = \
                (BCL2FASTQ_SEQNAMES["read1"], BCL2FASTQ_SEQNAMES["read2"],
                 BCL2FASTQ_SEQNAMES["sample_index"], BCL2FASTQ_SEQNAMES["barcode"])
            sample_names = read_chunk["sample_names"]
            sample_seq_bases = 0
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            for sample_name in sample_names:
                for seq_name in (r1_read, r2_read):
                    read_paths = find_func(path, seq_name, sample_name, lanes)
                    for read_fn in read_paths:
                        _, predicted_seq_bases = fastq_data_estimate(read_fn)
                        sample_seq_bases += predicted_seq_bases

            martian.log_info("Input data: Predict %f GB from %s" % (sample_seq_bases / 1e9, path))
            total_seq_bases += sample_seq_bases

            for sample_name in sample_names:
                r1_reads = find_func(path, r1_read, sample_name, lanes)
                r2_reads = find_func(path, r2_read, sample_name, lanes)

                # allow empty sample index case if all reads in lane are same sample
                sis = find_func(path, si_read, sample_name, lanes)
                if sis is None or len(sis) == 0:
                    sis = [None] * len(r1_reads)

                barcodes = find_func(path, bc_read, sample_name, lanes)
                if unbarcoded or len(barcodes) == 0:
                    martian.log_info('No barcodes available: ignoring sc processing')
                    barcodes = [None] * len(r1_reads)

                if not (len(r1_reads) == len(r2_reads) == len(barcodes)):
                    martian.log_info("Read 1 files: %s" % str(r1_reads))
                    martian.log_info("Read 2 files: %s" % str(r2_reads))
                    martian.log_info("Barcode files: %s" % str(barcodes))
                    martian.exit("Read1, Read2, and Barcode files are mismatched. Exiting pipeline")

                # calculate chunks
                for r1, r2, b, si in zip(r1_reads, r2_reads, barcodes, sis):
                    (flowcell, lane) = get_run_data(r1)
                    if sample_id is not None:
                        rg_string = ':'.join(str(item) for item in [sample_id, library_id, gem_group, flowcell, lane])
                    else:
                        rg_string = 'None:None:None:None:None'
                    new_chunk = {
                        'read1': r1, 'read2': r2, 'reads_interleaved': False, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

    martian.log_info("Input data: Predict %f total GB" % (total_seq_bases / 1e9))

    if len(chunks) == 0:
        martian.exit("No input FASTQs were found for the requested parameters.")

    if args.downsample is not None and args.downsample.get('subsample_rate', None) is None \
            and args.downsample.get('gigabases', None) is not None:
        global_subsample_rate = min(1.0, args.downsample['gigabases'] * 1e9 / total_seq_bases)
        martian.log_info("Input data downsampling: Requested: %.2f GB, Estimated Input: %.2f GB, Downsample Rate: %.3f"
                         % (args.downsample['gigabases'], total_seq_bases / 1e9, global_subsample_rate))
        for chunk in chunks:
            chunk['subsample_rate'] *= global_subsample_rate

    martian.log_info("Input reads: %s" % str(chunks))
    outs.chunks = chunks
    outs.read_groups = [rg for rg in read_groups]

    downsample_info = get_downsample_info(args.downsample, total_seq_bases)
    with open(outs.downsample_info, 'w') as downsample_out:
        tenkit.safe_json.dump_numpy(downsample_info, downsample_out)

    check_fastqs(outs.chunks)


def validate_input(args):
    """Does various parsing and checking of input arguments before we enter the main flow path
    """
    ok, msg = tk_preflight.check_gem_groups(args.sample_def)
    if not ok:
        martian.exit(msg)

    def check_key(n, dict_in, name, tys):
        if not name in dict_in:
            martian.exit("Entry %d in sample_def missing required field: %s" % (n, name))
        if not (type(dict_in[name]) in tys):
            martian.exit("Entry %d in sample_def for '%s' has incorrect type -- expecting %s, got %s" % (
                n, name, str(tys), type(dict_in[name])))

    for (idx, sample_item) in enumerate(args.sample_def):
        check_key(idx, sample_item, "read_path", [str, unicode])
        check_key(idx, sample_item, "lanes",  [list, type(None)])
        check_key(idx, sample_item, "gem_group", [int, type(None)])
        if args.input_mode == "BCL_PROCESSOR":
            check_key(idx, sample_item, "sample_indices", [list, type(None)])
        elif args.input_mode == "ILMN_BCL2FASTQ":
            check_key(idx, sample_item, "sample_names", [list, type(None)])

    if args.input_mode not in ["BCL_PROCESSOR", "ILMN_BCL2FASTQ"]:
        martian.throw("Unrecognized input_mode: %s" % args.input_mode)

    if args.downsample is not None:
        assert("gigabases" in args.downsample or "subsample_rate" in args.downsample)
        assert(not("gigabases" in args.downsample and "subsample_rate" in args.downsample))
        if 'subsample_rate' in args.downsample and args.downsample['subsample_rate'] is not None:
            assert(args.downsample['subsample_rate'] <= 1.0)


def get_downsample_info(downsample, total_seq_bases):
    """Returns information about input versus requested GB when downsampling
    """
    available_gb = total_seq_bases / 1e9
    requested_gb = None
    requested_rate = None
    downsample_succeeded = True

    if downsample is None:
        post_downsample_gb = available_gb
    else:
        if downsample.get('gigabases', None) is not None:
            requested_gb = float(downsample['gigabases'])
            post_downsample_gb = min(available_gb, requested_gb)
            if available_gb < requested_gb:
                martian.log_info("Downsample requested more GB than was available; will not downsample.")
                downsample_succeeded = False
        elif downsample.get('subsample_rate', None) is not None:
            requested_rate = float(downsample['subsample_rate'])
            post_downsample_gb = available_gb * requested_rate

    return {'available_gb': available_gb,
            'requested_gb': requested_gb,
            'requested_rate': requested_rate,
            'post_downsample_gb': post_downsample_gb,
            'downsample_succeeded': downsample_succeeded}


def check_fastqs(chunks):
    keys = ['read1', 'read2', 'barcode', 'sample_index']

    def check_fastq(fastq):
        # Check if fastq is readable
        if not os.access(fastq, os.R_OK):
            martian.exit("Do not have file read permission for FASTQ file: %s" % fastq)

        # Check if fastq is gzipped
        gzip_suffix = '.gz'
        is_gzip_fastq = True
        try:
            with gzip.open(fastq) as f:
                f.read(1)
        except:
            is_gzip_fastq = False

        if is_gzip_fastq and not fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is gzipped but filename does not have %s suffix: %s" % (fastq, gzip_suffix))
        if not is_gzip_fastq and fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is not gzipped but filename has %s suffix: %s" % (fastq, gzip_suffix))

    for chunk in chunks:
        for key in keys:
            fastq = chunk.get(key)
            if fastq is not None:
                check_fastq(fastq)


def fastq_data_estimate(fn, num_reads = 1000000):
    # Open reader
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
        is_gz = True
    else:
        reader = open(fn, 'r')
        is_gz = False

    gen = tk_fasta.read_generator_fastq(reader)
    rds = itertools.islice(gen, num_reads)

    input_lens = [(len(header) + len(r) + len(qual) + 4, len(r)) for (header, r, qual) in rds]
    total_seq_len = sum(x[1] for x in input_lens)
    total_data_len = sum(x[0] for x in input_lens)
    file_sz = os.path.getsize(fn)

    # NOTE: do not try and use the gzip footer containing the length of the compressed data
    # that only reflects the length of the final gzip block. A valid gzip file may have
    # many blocks, so that field cannot be relied upon.

    if is_gz:
        compressed_sz = reader.myfileobj.tell()
        predicted_sz = total_data_len / compressed_sz * file_sz
    else:
        predicted_sz = file_sz

    read_yield = len(input_lens) / total_data_len
    seq_yield = total_seq_len / total_data_len
    predicted_reads = read_yield * predicted_sz
    predicted_seq = seq_yield * predicted_sz

    # Log estimate of downsampling
    martian.log_info("Estimates for: %s" % fn)
    dbg_str =  "compressed_size: %.2f, predicted_size: %.2f" % \
               (file_sz / 1e9, predicted_sz / 1e9)
    martian.log_info(dbg_str)

    return (predicted_reads, predicted_seq)


def get_run_data(fn):
    """ Parse flowcell + lane from the first FASTQ record.
    NOTE: we don't check whether there are multiple FC / lanes in this file.
    """
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
    else:
        reader = open(fn, 'r')

    gen = tk_fasta.read_generator_fastq(reader)

    try:
        (name, seq, qual) = gen.next()
        try:
            flowcell, lane = re.split(':', name)[2:4]
        except ValueError:
            flowcell, lane = None, None
        return flowcell, lane
    except StopIteration:
        martian.exit("FASTQ is empty: %s" % fn)
