"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Mark PCR duplicates in a BAM file
"""
from __future__ import division

import itertools
import json
import math
import os
import numpy as np
import martian
import pickle
import pysam
import tenkit.bam as tk_bam
import tenkit.lane as tk_lane
from barcodes import get_read_barcode, load_barcode_whitelist, whitelist_mem_gb
from tools import ReferenceManager, create_bam_infile
from tools.io import index_bam, hierarchical_merge_bam, sort_bam, sort_bed, merge_keyed_bed
from tools.peaks import adjusted_position_pairs, is_chimeric_fragment
from collections import namedtuple, Counter
from constants import (SELF_FIVE_PRIME_POS_TAG, MATE_FIVE_PRIME_POS_TAG, LOW_MAPQ_THRESHOLD,
                       MATE_MAPPING_QUALITY_TAG, NO_BARCODE, TENX_PRODUCT_NAME)

__MRO__ = """
stage MARK_DUPLICATES(
    in  bam        input,
    in  string     reference_path,
    in  json       raw_barcode_counts,
    in  string     barcode_whitelist,
    out bam        output,
    out bam.bai    index,
    out csv        singlecell_mapping,
    out tsv.gz     fragments,
    out tsv.gz.tbi fragments_index,
    src py         "stages/processing/mark_duplicates",
) split (
    in  map        lane_map,
    in  string     chunk_start,
    in  string     chunk_end,
    in  int        chunk_num,
)
"""

# For optical duplicate detection
OPTICAL_DUPLICATE_DISTANCE = 100

# For diffusion duplicate detection, max distance over which diffusion is expected
MAX_DIFFUSION_DUP_DISTANCE = 25e3

SINGLE_CELL_KEYS = ["total", "duplicate", "chimeric", "unmapped",
                    "lowmapq", "mitochondrial", "passed_filters"]

ReadFootprint = namedtuple("ReadFootprint", ["barcode", "ref_id", "position", "mate_ref_id", "mate_position"])
NamedRead = namedtuple("NamedRead", ["footprint", "barcode", "read"])


# Generator utilities -- move to tenkit?
def consumer(func):
    """ decorator for initializing a generator consumer function """

    def start(*args, **kwargs):
        c = func(*args, **kwargs)
        c.next()
        return c

    return start


def broadcast(source, consumers):
    """ send each item in the source generator to each consumer in the list """
    for item in source:
        for c in consumers:
            c.send(item)

    for c in consumers:
        c.close()


def chunk_bound_func(read):
    if not read.is_unmapped:
        return read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG)
    else:
        return None


class DupSummary:
    def __init__(self, split_bcs, lane_coordinate_system, output_bam, output_tsv,
                 ref, bam_refs, priors=None, write_to_stdout=False):
        """ Summarize dups at a given subsampling rate, and barcode
            splitting policy.  If an open output_bam pysam.Samfile
            is passed, dups will be marked and reads will be written
            to output_bam """
        self.split_bcs = split_bcs
        if output_bam is None or output_tsv is None:
            raise ValueError("Must provide valid output paths")

        self.output_bam = output_bam
        self.output_tsv = output_tsv
        self.primary_contigs = ref.primary_contigs(allow_sex_chromosomes=True)
        self.mito_contigs = ref.non_nuclear_contigs()
        self.contig_lookup = bam_refs
        self.write_to_stdout = write_to_stdout
        self.bc_counts = {}
        self.contig_lengths = ref.get_contig_lengths()

        self.lane_coordinate_system = lane_coordinate_system

        # This is the raw counts of reads per barcode, used to discriminate when assigning a dupgroup to a barcode
        self.priors = priors

    def count_dups_by_distance(self, namedreads):
        """Count number of nearby duplicates in a set of reads.  A pair is counted as 1"""
        # Get (flowcell, lane, surface, swath, tile, x, y) tuples for each read
        read_locs = []
        for (footprint, barcode, read) in namedreads:
            read_loc = tk_lane.extract_read_position(read)
            if read_loc is not None:
                read_locs.append((read_loc, read))

        # Sort by flowcell_lane
        def flowcell_lane(read_loc):
            return "%s_%s" % (read_loc[0].flowcell, read_loc[0].lane)

        read_locs.sort(key=flowcell_lane)
        lane_groups = itertools.groupby(read_locs, flowcell_lane)

        opt_dups_found = 0  # really close dupes
        diff_dups_found = 0  # somewhat close dupes

        # Measure distances between all pairs in a lane
        for (lane, lane_reads) in lane_groups:
            lane_reads = list(lane_reads)

            layout = self.lane_coordinate_system.get_layout_for_read_loc(lane_reads[0][0])
            test_dups = layout.has_diffusion_duplicates(MAX_DIFFUSION_DUP_DISTANCE)

            if len(lane_reads) > 100:
                martian.log_info("Got dup cluster of size: %d" % len(lane_reads))
                first_read = lane_reads[0][1]
                martian.log_info("tid: %d, pos: %d, mapq: %d, seq: %s" % (
                    first_read.reference_id, first_read.reference_start, first_read.mapping_quality,
                    first_read.query_sequence))

            opt_dups = set()
            diff_dups = set()
            dump = []
            cmp_reads = min(200, len(lane_reads))
            lane_loc_coords = [self.lane_coordinate_system.convert_to_lane_coords(loc) for (loc, _) in lane_reads]
            for i in range(cmp_reads):
                loc1, read1 = lane_reads[i]
                lane_loc1 = lane_loc_coords[i]

                for j in range(i + 1, len(lane_reads)):
                    loc2, read2 = lane_reads[j]
                    lane_loc2 = lane_loc_coords[j]

                    dist = math.sqrt((lane_loc1[0] - lane_loc2[0]) ** 2 + (lane_loc1[1] - lane_loc2[1]) ** 2)
                    if test_dups and dist < MAX_DIFFUSION_DUP_DISTANCE:
                        diff_dups.add(j)
                        if self.write_to_stdout and j not in diff_dups:
                            dump.append(("%d\t" + ("%d\t" * 14)) % (dist,
                                                                    loc1.surface, loc1.swath, loc1.tile, loc1.x, loc1.y,
                                                                    lane_loc1[0], lane_loc1[1],
                                                                    loc2.surface, loc2.swath, loc2.tile, loc2.x, loc2.y,
                                                                    lane_loc2[0], lane_loc2[1]))

                    if dist < OPTICAL_DUPLICATE_DISTANCE:
                        opt_dups.add(j)

            if self.write_to_stdout and len(diff_dups) >= 2:
                for x in dump:
                    print ("%d\t%s" % (len(diff_dups), x))

            diff_dups_found += len(diff_dups)
            opt_dups_found += len(opt_dups)

        return opt_dups_found, diff_dups_found

    def process_read_block(self, reads):
        """dedups a block of reads, then writes them to output BAM in original order """
        read_tuples = []
        for read in reads:
            barcode = get_read_barcode(read)
            if barcode not in self.bc_counts:
                self.bc_counts[barcode] = {key: 0 for key in SINGLE_CELL_KEYS}

            if read.is_secondary:
                continue

            if read.is_unmapped or read.mate_is_unmapped or read.get_tag(SELF_FIVE_PRIME_POS_TAG) == read.get_tag(MATE_FIVE_PRIME_POS_TAG):
                # For unmapped pairs, key off of R1, and only report stats on R1
                read_key = read.is_read1
            else:
                # For mapped pairs, key so that the 5' most mate is primary
                read_key = read.get_tag(SELF_FIVE_PRIME_POS_TAG) < read.get_tag(MATE_FIVE_PRIME_POS_TAG)

            # We only need to dedup mapped pairs, but we output stats based on the read_key for consistency
            if read_key:
                self.bc_counts[barcode]['total'] += 1
                if read.is_unmapped or read.mate_is_unmapped:
                    self.bc_counts[barcode]['unmapped'] += 1

            if read.is_unmapped or read.mate_is_unmapped:
                continue

            # The footprint is what we form duplicate groups out of:  read barcode if split_bcs is set, and
            # read and mate contig IDs and 5' positions (as given by previously annotated tags)
            footprint = ReadFootprint(barcode if self.split_bcs else None,
                                      self.contig_lookup[read.reference_id], read.get_tag(SELF_FIVE_PRIME_POS_TAG),
                                      self.contig_lookup[read.next_reference_id], read.get_tag(MATE_FIVE_PRIME_POS_TAG))

            read_tuples.append(NamedRead(footprint, barcode, read))

        # Sort and then group by the read footprint.  Note that the sort is necessary to group all reads with the
        # same footprint.
        read_tuples.sort(key=lambda x: x.footprint)
        dup_groups = itertools.groupby(read_tuples, lambda x: x.footprint)

        for (footprint, dup_group) in dup_groups:
            dup_group = list(dup_group)
            total_dups = len(dup_group)
            contig = footprint.ref_id
            if total_dups > 1:
                optical_dups, diffusion_dups = self.count_dups_by_distance(dup_group)
            else:
                optical_dups = 0
                diffusion_dups = 0
            non_proximal_dups = total_dups - max(diffusion_dups, optical_dups)

            dup_group_barcodes = Counter()
            for namedread in dup_group[:non_proximal_dups]:
                if namedread.barcode is not None:
                    dup_group_barcodes[namedread.barcode] += 1
            if not dup_group_barcodes:
                most_common_barcode = None
            else:
                # Use raw read counts per barcode to break ties in determining the best barcode
                if self.priors is None:
                    most_common_barcode = dup_group_barcodes.most_common(1)[0][0]
                else:
                    max_count = max(dup_group_barcodes.values())
                    common_barcodes = [bc for bc, count in dup_group_barcodes.iteritems()
                                       if count == max_count]
                    most_common_barcode = max(common_barcodes, key=lambda bc: self.priors[bc])

            # Identify the unique duplicate out of the group as the one with the minimum query name to be consistent
            # between this dup group and its read pairs
            unique_dup_index = min([(i, dup_group[i].read.query_name) for i in range(total_dups)
                                    if dup_group[i].barcode == most_common_barcode], key=lambda x: x[1])[0]

            unique_read = dup_group[unique_dup_index].read

            # NOTE: this means that number of dups as per BAM tag is slightly different from number of dups counted
            # by discarding lowmapq, mito or chimeric fragments.
            for i in range(0, total_dups):
                dup_group[i].read.is_duplicate = True
            unique_read.is_duplicate = False

            fragment_mapq = min(unique_read.mapping_quality, unique_read.get_tag(MATE_MAPPING_QUALITY_TAG))
            primary_contigs_set = set(self.primary_contigs)
            # Make sure we only output data for one read dup group of each read pair
            equal_positions = footprint.position == footprint.mate_position
            if footprint.position <= footprint.mate_position:
                if fragment_mapq <= LOW_MAPQ_THRESHOLD:
                    # Count chimerically mapped fragments with identical 5' tags for each read only once
                    for read in (dup_group[i] for i in range(total_dups) if not equal_positions or dup_group[i].read.is_read1):
                        self.bc_counts[read.barcode]["lowmapq"] += 1
                elif contig in self.mito_contigs:
                    for read in (dup_group[i] for i in range(total_dups) if not equal_positions or dup_group[i].read.is_read1):
                        self.bc_counts[read.barcode]["mitochondrial"] += 1
                elif is_chimeric_fragment(unique_read) or contig not in primary_contigs_set:
                    # Note that we've added fragments mapping to supplementary contigs here
                    for read in (dup_group[i] for i in range(total_dups) if not equal_positions or dup_group[i].read.is_read1):
                        self.bc_counts[read.barcode]["chimeric"] += 1
                else:
                    self.bc_counts[most_common_barcode]["passed_filters"] += 1
                    for dup_read in (dup_group[i] for i in range(total_dups) if i != unique_dup_index):
                        self.bc_counts[dup_read.barcode]["duplicate"] += 1
                    if most_common_barcode is not None:
                        # Write out the output fragments
                        start, stop = adjusted_position_pairs(unique_read)
                        if start is not None and stop is not None:
                            start = max(0, start)
                            stop = min(stop, self.contig_lengths[contig])
                            self.output_tsv.write('{contig}\t{start}\t{stop}\t'
                                                  '{most_common_barcode}\t'
                                                  '{non_proximal_dups}\n'.format(**locals()))

        for read in reads:
            self.output_bam.write(read)

    @consumer
    def read_consumer(self):
        # bam is sorted by SELF_FIVE_PRIME_POS tag, chrom and pos.
        current_bam_key = (-1, -1)
        current_reads = []
        try:
            while True:
                # accept the next read
                read = (yield)

                new_bam_key = (read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG))

                # If the dup group gets extremely large we can run out of memory.
                # Process things in groups of 500K to prevent memory blow-up
                # May cause us to miss a few dups, but it doesn't really matter in these crazy regions
                if new_bam_key != current_bam_key or len(current_reads) > 500000:
                    process_reads = current_reads
                    current_reads = []
                    current_bam_key = new_bam_key

                    if len(process_reads) > 0:
                        self.process_read_block(process_reads)

                # accumulate block of reads with same start position
                current_reads.append(read)

        except GeneratorExit:
            # Finish up final batch
            self.process_read_block(current_reads)
            return


def split(args):
    # Chunk bam to get 1GB per chunk
    bam_in = create_bam_infile(args.input)
    bam_chunk_size_disk = 0.75
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_func, chunk_size_gb=bam_chunk_size_disk)

    for chunk in chunk_defs:
        chunk['__mem_gb'] = 4
        chunk['__vmem_gb'] = 5 + int(np.ceil(2 * whitelist_mem_gb(args.barcode_whitelist) + bam_chunk_size_disk * 10))

    lane_coord_sys = tk_lane.LaneCoordinateSystem()

    # Reopen BAM for estimating tile extents
    bam_in = create_bam_infile(args.input)
    lane_coord_sys.estimate_tile_extents(bam_in)
    for cnum, chunk in enumerate(chunk_defs):
        chunk['lane_map'] = lane_coord_sys.to_dict()
        chunk['chunk_num'] = cnum

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 8, '__threads': 4}}


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    # Merge the output bam files with duplicates marked
    hierarchical_merge_bam([c.output for c in chunk_outs], outs.output, tag=None,
                           threads=martian.get_threads_allocation())
    outs.index = index_bam(outs.output, martian.get_threads_allocation())

    # Merge the barcode counts from each chunk and write out our singlecell_mapping file
    barcode_whitelist = load_barcode_whitelist(args.barcode_whitelist, ordered=True)
    sorted_barcodes = []
    if args.raw_barcode_counts is not None:
        with open(args.raw_barcode_counts, 'r') as infile:
            raw_counts = json.load(infile)
        sorted_barcodes = ['{}-{}'.format(barcode, gem_group)
                           for gem_group in raw_counts
                           for barcode in sorted(barcode_whitelist)]
    barcode_counts = {}
    for chunk in chunk_outs:
        with open(chunk.singlecell_mapping, 'r') as infile:
            chunk_counts = pickle.load(infile)
        for barcode, count_dict in chunk_counts.iteritems():
            if barcode not in barcode_counts:
                barcode_counts[barcode] = Counter()
            barcode_counts[barcode] += Counter(count_dict)

    with open(outs.singlecell_mapping, 'w') as outfile:
        outfile.write("barcode,")
        outfile.write(",".join(SINGLE_CELL_KEYS))
        outfile.write("\n")
        if None in barcode_counts:
            outfile.write("{},".format(NO_BARCODE))
            outfile.write(",".join([str(barcode_counts[None][key]) for key in SINGLE_CELL_KEYS]))
            outfile.write("\n")
        for barcode in (bc for bc in sorted_barcodes if bc in barcode_counts):
            outfile.write("{},".format(barcode))
            outfile.write(",".join([str(barcode_counts[barcode][key]) for key in SINGLE_CELL_KEYS]))
            outfile.write("\n")

    # Merge the fragment file
    base_file, extension = os.path.splitext(outs.fragments)
    if not extension == '.gz':
        raise ValueError('Expecting compressed file output')
    input_tsvs = [str(chunk.fragments) for chunk in chunk_outs]
    merge_keyed_bed(input_tsvs, base_file, threads=martian.get_threads_allocation())
    if os.path.getsize(base_file) == 0:
        outs.fragments = None
        outs.fragments_index = None
        return

    # N.B. tabix_index will automatically compress the input file, adding the .gz suffix
    pysam.tabix_index(base_file, preset='bed', index=outs.fragments_index)


def main(args, outs):
    """Mark exact duplicate reads in the output BAM file while also writing out some summary statistics.
    PCR duplicates have the same read1 start site and read2 start site.
    """
    args.coerce_strings()
    outs.coerce_strings()

    # Chunk output doesn't get indexed
    outs.fragments_index = None
    outs.index = None

    # Pull in prior likelihoods for barcodes
    raw_barcode_abundance = None
    barcode_whitelist = load_barcode_whitelist(args.barcode_whitelist)
    if args.raw_barcode_counts is not None and barcode_whitelist is not None:
        with open(args.raw_barcode_counts, 'r') as infile:
            raw_counts = json.load(infile)
        raw_barcode_abundance = {'{}-{}'.format(barcode, gem_group): count
                                 for gem_group, subdict in raw_counts.iteritems()
                                 for barcode, count in zip(barcode_whitelist, subdict['bc_counts'])}

    bam_in = create_bam_infile(args.input)
    bam_refs = bam_in.references

    bam_prefix, ext = os.path.splitext(outs.output)
    raw_bam_file = martian.make_path(bam_prefix + '_five_prime_pos_sorted' + ext)

    frag_prefix, ext = os.path.splitext(outs.fragments)
    raw_frag_file = martian.make_path(frag_prefix + '_raw' + ext)

    # only write CO line for one chunk, so we don't have duplicates after samtools merge
    if args.chunk_num == 0:
        COs = ['10x_bam_to_fastq:R1(SEQ:QUAL,TR:TQ)',
               '10x_bam_to_fastq:R2(SEQ:QUAL,TR:TQ)',
               '10x_bam_to_fastq:I1(BC:QT)',
               '10x_bam_to_fastq:I2(CR:CY)',
               '10x_bam_to_fastq_seqnames:R1,R3,I1,R2']
    else:
        COs = None

    bam_out, _ = tk_bam.create_bam_outfile(raw_bam_file, None, None, template=bam_in,
                                           pgs=[tk_bam.make_pg_header(martian.get_pipelines_version(),
                                                                      "mark_duplicates",
                                                                      TENX_PRODUCT_NAME)],
                                           cos=COs)
    fragments_out = open(raw_frag_file, 'w')
    bam_in.reset()

    # Ensure the summary key indicates what kind of dup marking was actually performed.
    lane_coord_sys = tk_lane.LaneCoordinateSystem.from_dict(args.lane_map)
    reference_manager = ReferenceManager(args.reference_path)
    summarizer = DupSummary(split_bcs=False,
                            lane_coordinate_system=lane_coord_sys,
                            output_bam=bam_out,
                            output_tsv=fragments_out,
                            ref=reference_manager,
                            bam_refs=bam_refs,
                            priors=raw_barcode_abundance)

    # Now broadcast the selected reads to the summarizers
    consumers = [summarizer.read_consumer()]
    source = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    broadcast(source, consumers)

    # Close outfiles
    bam_out.close()
    fragments_out.close()

    # Feed the chunk barcode_counts data back to join()
    with open(outs.singlecell_mapping, 'w') as outfile:
        pickle.dump(summarizer.bc_counts, outfile)

    # Sort the output bam & tsv files
    sort_bam(raw_bam_file, outs.output, threads=martian.get_threads_allocation())
    sort_bed(raw_frag_file, outs.fragments, genome=reference_manager.fasta_index,
             threads=martian.get_threads_allocation(), leave_key=True)
