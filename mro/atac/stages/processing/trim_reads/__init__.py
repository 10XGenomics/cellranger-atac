"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Consume multiple fastq files and spit out a pooled fastq file with targeted adapters removed from the ends of the reads.
"""
from __future__ import division

import itertools
import json
import random
from subprocess import PIPE

from barcodes import load_barcode_whitelist
import martian
import numpy as np
import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json
import tenkit.stats as tk_stats
import tenkit.seq as tk_seq
from tenkit.constants import WHITELIST_TO_LOT_MAP
from tools import open_maybe_gzip
from cellranger.fastq import infer_barcode_reverse_complement

__MRO__ = """
stage TRIM_READS(
    in  map[]  chunks,
    in  string barcode_whitelist,
    in  int    max_read_num,
    in  map    trim_def,
    in  map    adapters,
    out map[]  chunks,
    out json   bc_counts,
    out json   lot_info,
    out json   read_counts,
    src py     "stages/processing/trim_reads",
) split using (
    in  map    chunk,
)
"""

MEMORY_IN_GB = 3
MAX_BARCODE_ERROR_RATE = 0.97


def identify_gelbead_lot(barcode_hist, lot_to_barcodes, min_frac=0.95, min_counts=1000):
    # flip from {lot: [bcs]} to {bc: lot} for easier lookup
    barcode_to_lot = {barcode: lot
                      for lot, barcode_list in lot_to_barcodes.items()
                      for barcode in barcode_list}
    lot_counts = {lot: 0 for lot in lot_to_barcodes}
    part_a_len = len(barcode_to_lot.keys()[0])  # infer length of part A

    # count the total barcode observations for each lot
    for barcode, count in barcode_hist.items():
        part_a = barcode[:part_a_len]
        lot = barcode_to_lot.get(part_a)
        if lot is not None:
            lot_counts[lot] += count

    # determine best match
    best_lot = max(lot_counts, key=lambda lot: lot_counts[lot])
    best_counts = lot_counts[best_lot]
    total_counts = sum(lot_counts.values())

    best_frac = best_counts / total_counts if total_counts > 0 else 0.0

    if best_frac >= min_frac and total_counts >= min_counts:
        result_lot = best_lot
        result_conf = "confident"
    elif total_counts < min_counts:
        result_lot = None
        result_conf = "insufficient data"
    else:
        result_lot = None
        result_conf = "ambiguous"

    return result_lot, result_conf, lot_counts


def run_cutadapt_single_end(in_reads_fn, out_reads_fn, trim_info_fn, trim_def, adapters, read_id="R1"):
    """Calls cutadapt in single-end mode using the settings in trim_def[read_id]
    """
    filter_output = trim_def["discard_untrimmed"]
    if "trim_length" in trim_def:
        fixed_length = trim_def["trim_length"]
    else:
        fixed_length = None

    martian.log_info("Trim definition provided:\n{}".format(trim_def))
    martian.log_info("(Using info for read {}".format(read_id))
    martian.log_info("Adapter sequences provided:\n{}".format(adapters))

    seqs_to_trim = {}
    for direction in ["5prime", "3prime"]:
        if read_id in trim_def and direction in trim_def[read_id]:
            seqs_to_trim[direction] = "".join(adapters[idx] for idx in trim_def[read_id][direction])
        else:
            seqs_to_trim[direction] = None

    cmd = ["cutadapt"]
    if fixed_length is not None:
        cmd.extend(["--length", "{}".format(fixed_length)])

    start_trim = seqs_to_trim["5prime"]
    end_trim = seqs_to_trim["3prime"]
    read_trim = start_trim is not None or end_trim is not None
    if start_trim is not None and end_trim is not None:
        # This is a linked adapter trim that anchors the 5prime end adapter to the beginning of the read
        # and lets the 3prime adapter float
        cmd.extend(["-a", "{}...{}".format(start_trim, end_trim)])
    elif start_trim is not None:
        # Just the anchored 5prime end adapter
        cmd.extend(["-g", "^{}".format(start_trim)])
    elif end_trim is not None:
        # Just the floating 3prime end adapter
        cmd.extend(["-a", end_trim])

    if filter_output and read_trim:
        cmd.append("--discard-untrimmed")

    cmd.extend(["--info-file", trim_info_fn])
    cmd.extend(["-o", out_reads_fn])
    cmd.append(in_reads_fn)

    martian.log_info("Cutadapt command: \n{}".format(" ".join(cmd)))

    process = tk_subproc.Popen(cmd, stdin=None, stdout=PIPE, stderr=PIPE)
    (stdout, stderr) = process.communicate()
    if process.returncode != 0:
        martian.log_info("Error while running cutadapt: \n{}".format(stderr))
        raise ValueError("Cutadapt failed")

    martian.log_info("Cutadapt output: \n{}".format(stdout))

    input_read_pairs, output_read_pairs = None, None
    for line in stdout.split("\n"):
        if line.startswith("Total reads processed:"):
            input_read_pairs = int(line.split(":")[1].replace(",", ""))
        if line.startswith("Reads written (passing filters):"):
            output_read_pairs = int(line.split(":")[1].split("(")[0].replace(",", ""))

    return input_read_pairs, output_read_pairs


def split(args):
    """Chunking is handled by a prior stage (SETUP_CHUNKS), so we just have to set them up.
    """
    chunks = [{'chunk': chunk, "__mem_gb": MEMORY_IN_GB}
              for chunk in args.chunks]
    return {'chunks': chunks}


def join(args, outs, chunk_defs, chunk_outs):
    final_chunks = []

    for cl in chunk_outs:
        final_chunks.extend(cl.chunks)

    outs.chunks = final_chunks
    valid_counts = [c.bc_counts for c in chunk_outs if c.bc_counts is not None]

    # No counts if there's no whitelist or actual counts
    if args.barcode_whitelist is None or len(valid_counts) == 0:
        outs.bc_counts = None
        outs.lot_info = None
        return

    bc_counts = {}
    read_counts = {}
    for (c_out, c_def) in zip(chunk_outs, chunk_defs):
        # Sum up total and trimmed read counts
        with open(c_out.read_counts) as f:
            r = json.load(f)
        for key in ['filtered_read_pairs', 'total_read_pairs']:
            read_counts[key] = read_counts.get(key, 0) + r[key]

        # Sum up barcode counts
        gem_group = c_def.chunk['gem_group']
        if c_out.bc_counts is not None:
            with open(c_out.bc_counts) as f:
                r = json.load(f)
            gg_result = bc_counts.setdefault(gem_group, {'bad_bc_count': 0, 'bc_counts': None})
            gg_result['bad_bc_count'] += r['bad_bc_count']
            if gg_result['bc_counts'] is None:
                gg_result['bc_counts'] = np.array(r['bc_counts'], dtype=np.int32)
            else:
                gg_result['bc_counts'] += np.array(r['bc_counts'], dtype=np.int32)

    total_counts = 0
    total_errors = 0
    for gg in bc_counts.keys():
        rgg = bc_counts[gg]
        rgg['bc_error_rate'] = tk_stats.robust_divide(float(rgg['bad_bc_count']),
                                                      float(rgg['bad_bc_count'] + rgg['bc_counts'].sum()))
        total_counts += float(rgg['bad_bc_count'] + rgg['bc_counts'].sum())
        total_errors += float(rgg['bad_bc_count'])

    # Hardcoded bail-out if the BC-error rate is extremely high
    bc_error_rate = total_errors / total_counts
    if bc_error_rate > MAX_BARCODE_ERROR_RATE:
        martian.exit("Extremely high rate of incorrect barcodes observed (%.2f %%). "
                     "Check that input is 10x Chromium data, "
                     "and that there are no missing cycles in first 16 bases of the index read I2." %
                     (bc_error_rate * 100.0))

    # possibly do lot detection
    lot_detection = {}
    lot_map = WHITELIST_TO_LOT_MAP.get(args.barcode_whitelist, None)
    if lot_map is not None:
        # get BC counts histogram
        # for now, just sum over all gem groups
        bc_seq = sorted(list(load_barcode_whitelist(args.barcode_whitelist)))
        bc_cts = np.sum([ggr['bc_counts'] for ggr in bc_counts.values()], axis=0)
        bc_hist = {seq: cts for seq, cts in zip(bc_seq, bc_cts)}

        (gelbead_lot, gelbead_lot_confidence, gelbead_lot_counts) = identify_gelbead_lot(bc_hist, lot_map)
        # only report on lots with nonzero counts
        gelbead_lot_counts_nonzero = {lot: count for lot, count in gelbead_lot_counts.items() if count > 0}

        lot_detection['gelbead_lot'] = gelbead_lot
        lot_detection['gelbead_lot_confidence'] = gelbead_lot_confidence
        lot_detection['gelbead_lot_counts'] = gelbead_lot_counts_nonzero

        martian.log_info("Gelbead lot detected: %s, reason (if None): %s" % (gelbead_lot, gelbead_lot_confidence))

    if outs.lot_info is not None:
        with open(outs.lot_info, 'w') as f:
            tenkit.safe_json.dump_numpy(lot_detection, f, pretty=True)

    if outs.bc_counts is not None:
        with open(outs.bc_counts, 'w') as f:
            tenkit.safe_json.dump_numpy(bc_counts, f)

    with open(outs.read_counts, 'w') as f:
        tenkit.safe_json.dump_numpy(read_counts, f, pretty=True)


def read_match(read, other):
    """Compares two read names, checking to see if they're part of the same read pair or index reads
    """
    return read[0].split()[0] == other[0].split()[0]

def get_read_name(read):
    """Get read names split out"""
    return read[0].split()[0]

def read_generator_trim_info(trim_info_file, paired_end=False):
    """Takes the output from cutadapt --info-file and gives an iterator over the trimmed sequences/qualities
    compatible with that produced by tenkit.fasta.read_generator_fastq.
    N.B.: assumes that sequences are trimmed from the 3' end only.

    Example cutadapt output:
    (no adapter match)
    A00228:277:HFKLHDMXX:1:1110:17716:1078 1:N:0:0	-1	CTTCTAGTCATAGCATCACTCTATTTTTTTTTTTTTTTTTTGAGACAGAG	F::FFF:FFFF:FFFFFFFFFFFF:FFF:FF:FFFFFFFF,::,FFFFF:
    (adapter matched to 16 bp on 3' end)
    A00228:277:HFKLHDMXX:1:1160:21070:30467 1:N:0:0	0	34	50	AGACAGAGTCTTGCTCTGTGGCCTAGGCTGGAGT	CTGTCTCTTATACACA		1	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	FFFFFFFFFFFFFFFF
    """
    # fewer than this many columns means no adapter match
    CUTADAPT_NCOLS_MATCH = 11

    # column numbers in cutadapt --info-file output:
    CUTADAPT_READ_NAME = 0
    CUTADAPT_SEQ_MATCHED = 5    # sequence matched to the adapter sequence
    CUTADAPT_SEQ_RIGHT = 6      # sequence to the right (3' direction) of adapter match
    CUTADAPT_QUAL_MATCHED = 9   # quality scores for sequence matched to adapter
    CUTADAPT_QUAL_RIGHT = 10    # quality scores to the right of adapter match

    for line1 in trim_info_file:
        if paired_end:
            try:
                line2 = trim_info_file.next()
            except StopIteration:
                raise EOFError("File output by `cutadapt --info-file` had an odd number of lines.")

        read1_info = line1.strip('\n').split('\t')
        read1_name = read1_info[CUTADAPT_READ_NAME]
        if len(read1_info) == CUTADAPT_NCOLS_MATCH:    # found adapter match
            read1_trim_seq = read1_info[CUTADAPT_SEQ_MATCHED] + read1_info[CUTADAPT_SEQ_RIGHT]
            read1_trim_qual = read1_info[CUTADAPT_QUAL_MATCHED] + read1_info[CUTADAPT_QUAL_RIGHT]
        else:                   # no adapter match
            read1_trim_seq = ''
            read1_trim_qual = ''

        if paired_end:
            read2_info = line2.strip('\n').split('\t')
            read2_name = read2_info[CUTADAPT_READ_NAME]
            if len(read2_info) == CUTADAPT_NCOLS_MATCH:    # found adapter match
                read2_trim_seq = read2_info[CUTADAPT_SEQ_MATCHED] + read2_info[CUTADAPT_SEQ_RIGHT]
                read2_trim_qual = read2_info[CUTADAPT_QUAL_MATCHED] + read2_info[CUTADAPT_QUAL_RIGHT]
            else:                   # no adapter match
                read2_trim_seq = ''
                read2_trim_qual = ''
            yield (read1_name, read1_trim_seq, read1_trim_qual,
                   read2_name, read2_trim_seq, read2_trim_qual)
        else:
            yield (read1_name, read1_trim_seq, read1_trim_qual)


def main(args, outs):
    """ Trim the reads in a series of fasta files """
    chunk = args.chunk
    subsample_rate = chunk['subsample_rate']
    have_barcode = chunk['barcode'] is not None
    have_sample_index = chunk['sample_index'] is not None

    # STEP 1:  We run the R1/R2 reads through cutadapt, writing them to a temporary file with appropriate adapters
    # trimmed, optionally filtering out reads where adapters weren't found
    interleaved = chunk['read2'] is None
    # can't do discard_untrimmed because we're running cutadapt in single-end mode
    if args.trim_def['discard_untrimmed']:
        martian.exit("discard_untrimmed was set in trim_def")
    if interleaved:
        trimmed_reads = martian.make_path("trimmed_reads.fastq.gz")
        trim_info_fn = martian.make_path("trim_info.txt.gz")
        initial_read_pairs, trimmed_read_pairs = run_cutadapt_single_end(chunk['read1'], trimmed_reads,
                                                                         trim_info_fn, args.trim_def, args.adapters)
    else:
        trimmed_r1 = martian.make_path("trimmed_r1.fastq.gz")
        trimmed_r2 = martian.make_path("trimmed_r2.fastq.gz")
        trim_info_r1_fn = martian.make_path("trim_info_r1.txt.gz")
        trim_info_r2_fn = martian.make_path("trim_info_r2.txt.gz")
        initial1, trimmed1 = run_cutadapt_single_end(chunk['read1'], trimmed_r1,
                                                     trim_info_r1_fn, args.trim_def, args.adapters, read_id="R1")
        initial2, trimmed2 = run_cutadapt_single_end(chunk['read2'], trimmed_r2,
                                                     trim_info_r2_fn, args.trim_def, args.adapters, read_id="R2")
        initial_read_pairs = initial1 + initial2
        trimmed_read_pairs = trimmed1 + trimmed2
        if initial1 != initial2:
            martian.exit("Input fastq files for R1 and R2 are not the same length")
        if trimmed1 != trimmed2:
            raise ValueError("Cutadapt produced differing numbers of reads for R1 and R2")

    # STEP 2:  We run through the trimmed R1/R2 reads along with sample index and barcode reads, chunking into files of
    # max_read_num reads or less, and skipping sample index/barcode reads that don't match the trimmed & filtered R1/R2
    # reads
    max_read_num = args.max_read_num
    file_number = 1

    # open the available input read files and get the iterator over them
    if interleaved:
        reads_in = open_maybe_gzip(trimmed_reads, 'r')
        read_iter = tk_fasta.read_generator_fastq(reads_in, paired_end=True)
        trim_info = open_maybe_gzip(trim_info_fn, 'r')
        trim_iter = read_generator_trim_info(trim_info, paired_end=True)
    else:
        r1_in = open_maybe_gzip(trimmed_r1, 'r')
        r2_in = open_maybe_gzip(trimmed_r2, 'r')
        read_iter = ((r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])
                     for r1, r2 in itertools.izip_longest(tk_fasta.read_generator_fastq(r1_in),
                                                          tk_fasta.read_generator_fastq(r2_in)))
        trim_info_r1 = open_maybe_gzip(trim_info_r1_fn, 'r')
        trim_info_r2 = open_maybe_gzip(trim_info_r2_fn, 'r')
        trim_iter = (t1 + t2 for t1, t2 in itertools.izip(read_generator_trim_info(trim_info_r1),
                     read_generator_trim_info(trim_info_r2)))

    # open output read file, which will be interleaved
    read_name = martian.make_path("read{}.fastq.lz4".format(file_number))
    out_readfiles = [read_name]
    out_read_fastq = open_maybe_gzip(read_name, 'w')

    # open trimmed read file, which will be interleaved
    trim_out_name = martian.make_path("TRIM{}.fastq.lz4".format(file_number))
    out_trimfiles = [trim_out_name]
    out_trim_fastq = open_maybe_gzip(trim_out_name, 'w')

    if args.barcode_whitelist is None:
        outs.bc_counts = None
        barcode_indices = None
    else:
        barcode_whitelist = sorted(list(load_barcode_whitelist(args.barcode_whitelist)))
        barcode_indices = {bc: idx for (idx, bc) in enumerate(barcode_whitelist)}
        bc_counts = np.zeros(len(barcode_whitelist), dtype=np.int32)
        bad_count = 0

    # open barcode file if there is one
    if have_barcode:
        bc_name = martian.make_path("BC{}.fastq.lz4".format(file_number))
        out_bc_fastq = open_maybe_gzip(bc_name, 'w')
        out_barcodefiles = [bc_name]
        barcode_read = None
        bc_in = open_maybe_gzip(chunk['barcode'], 'r')
        bc_iter = tk_fasta.read_generator_fastq(bc_in)
        # Determine if barcode sequences need to be reverse complemented.
        with open_maybe_gzip(chunk['barcode'], 'r') as bc_in2:
            bc_iter2 = tk_fasta.read_generator_fastq(bc_in2)
            barcode_whitelist = load_barcode_whitelist(args.barcode_whitelist)
            barcode_rc = infer_barcode_reverse_complement(barcode_whitelist, bc_iter2)
    else:
        out_barcodefiles = [None]
        outs.bc_counts = None

    # open sample_index file if there is one
    if have_sample_index:
        si_name = martian.make_path("SI{}.fastq.lz4".format(file_number))
        out_si_fastq = open_maybe_gzip(si_name, 'w')
        si_in = open_maybe_gzip(chunk['sample_index'], 'r')
        sample_index_read = None
        si_iter = tk_fasta.read_generator_fastq(si_in)
        out_sampleindex_files = [si_name]
    else:
        out_sampleindex_files = [None]

    read_num = 0
    random.seed(0)
    for (read, trim) in itertools.izip(read_iter, trim_iter):
        # Downsample (other than the first read).  Note we've set a fixed seed to make this deterministic.
        if read_num > 0 and random.random() > subsample_rate:
            continue

        # Now we need to step through the barcode and sample index reads to find the matching reads
        if have_barcode:
            try:
                while barcode_read is None or not read_match(read, barcode_read):
                    barcode_read = bc_iter.next()
                # reverse complement if all barcodes are RC-ed
                if barcode_rc:
                    barcode_read = (barcode_read[0], tk_seq.get_rev_comp(barcode_read[1]), barcode_read[2][::-1])
            except StopIteration:
                raise ValueError("Couldn't find barcode read matching {}".format(get_read_name(read)))
        if have_sample_index:
            try:
                while sample_index_read is None or not read_match(read, sample_index_read):
                    sample_index_read = si_iter.next()
            except StopIteration:
                raise ValueError("Couldn't find sample index read matching {}".format(get_read_name(read)))

        (name1, seq1, qual1, name2, seq2, qual2) = read
        (tr_name1, tr_seq1, tr_qual1, tr_name2, tr_seq2, tr_qual2) = trim

        read_num += 1
        if read_num > max_read_num:
            read_num = 1
            file_number += 1
            read_name = martian.make_path("read{}.fastq.lz4".format(file_number))
            out_read_fastq.close()
            out_read_fastq = open_maybe_gzip(read_name, 'w')
            out_readfiles.append(read_name)

            trim_out_name = martian.make_path("TRIM{}.fastq.lz4".format(file_number))
            out_trim_fastq.close()
            out_trim_fastq = open_maybe_gzip(trim_out_name, 'w')
            out_trimfiles.append(trim_out_name)

            if have_barcode:
                bc_name = martian.make_path("BC{}.fastq.lz4".format(file_number))
                out_bc_fastq.close()
                out_bc_fastq = open_maybe_gzip(bc_name, 'w')
                out_barcodefiles.append(bc_name)
            else:
                out_barcodefiles.append(None)

            if have_sample_index:
                si_name = martian.make_path("SI{}.fastq.lz4".format(file_number))
                out_si_fastq.close()
                out_si_fastq = open_maybe_gzip(si_name, 'w')
                out_sampleindex_files.append(si_name)
            else:
                out_sampleindex_files.append(None)

        if have_barcode:
            barcode_seq = barcode_read[1]
            barcode_qual = barcode_read[2]
            if barcode_indices is not None:
                idx = barcode_indices.get(barcode_seq)
                if idx is not None:
                    bc_counts[idx] += 1
                else:
                    bad_count += 1

            tk_fasta.write_read_fastq(out_bc_fastq, barcode_read[0], barcode_seq, barcode_qual)
        if have_sample_index:
            tk_fasta.write_read_fastq(out_si_fastq, sample_index_read[0], sample_index_read[1], sample_index_read[2])

        tk_fasta.write_read_fastq(out_read_fastq, name1, seq1, qual1)
        tk_fasta.write_read_fastq(out_read_fastq, name2, seq2, qual2)

        tk_fasta.write_read_fastq(out_trim_fastq, tr_name1, tr_seq1, tr_qual1)
        tk_fasta.write_read_fastq(out_trim_fastq, tr_name2, tr_seq2, tr_qual2)

    if interleaved:
        reads_in.close()
    else:
        r1_in.close()
        r2_in.close()

    if have_barcode:
        out_bc_fastq.close()
        # Only emit BC counts if we had a whitelist
        if outs.bc_counts is not None:
            result = {}
            result['bad_bc_count'] = bad_count
            result['bc_counts'] = list(bc_counts)
            with open(outs.bc_counts, 'w') as bc_counts_out:
                tenkit.safe_json.dump_numpy(result, bc_counts_out)

    with open(outs.read_counts, 'w') as outfile:
        read_counts = {'total_read_pairs': initial_read_pairs,
                       'filtered_read_pairs': trimmed_read_pairs}
        tenkit.safe_json.dump_numpy(read_counts, outfile)

    if have_sample_index:
        out_si_fastq.close()
    out_read_fastq.close()
    out_trim_fastq.close()

    outs.chunks = [{'read1': r,  # output chunked trimmed read file
                    'read2': None,
                    'trim': t,          # output chunked trim file
                    'barcode': bc,  # output chunked barcode file
                    'sample_index': si,  # output chunked sample index file
                    'barcode_reverse_complement': False,  # we always keep BC in correct orientation
                    'reads_interleaved': True,
                    'gem_group': chunk['gem_group'],
                    'read_group': chunk['read_group']}
                   for (r, t, bc, si) in
                   zip(out_readfiles, out_trimfiles, out_barcodefiles, out_sampleindex_files)]
