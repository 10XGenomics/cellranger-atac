import os
import os.path
import resource
import shutil
import subprocess
import sys
import pysam
import pandas as pd
import numpy as np
import tenkit.log_subprocess as tk_subproc
import tools.regions as regtools

from collections import namedtuple, Counter
from cellranger.h5_constants import GZIP_SUFFIX, LZ4_SUFFIX
from tools import ReferenceManager

def parse_aggr_csv(infile, whitelist=["library_id", "fragments", "cells"], blacklist=["batch"]):
    aggr_df = pd.read_csv(infile, sep=',')
    nchunks = len(aggr_df)
    library_info = {}
    if aggr_df.isnull().values.any():
        return None, None, "ill-formatted file {}, with either empty strings as values or missing values in general".format(infile)

    aggr_keys = aggr_df.columns.values
    libraries = set()
    if len(set(whitelist).intersection(set(aggr_keys))) == 0:
        return None, None, "Did you forget to include a header line in {}?".format(infile)

    whitespace = False
    for key in aggr_keys:
        whitespace |= (key.lstrip().rstrip() != key)

    for key in whitelist:
        if key not in aggr_keys:
            msg = "Your header row is missing a required field: {}.".format(key)
            if whitespace:
                return None, None, msg + "At least one of the header fields contains whitespace which may be unintended. Whitespaces will be interpreted as part of input"
            return None, None, msg
    for key in aggr_keys:
        if blacklist is not None and key in blacklist:
            return None, None, "Invalid/unsupported header field: {}".format(key)
        for num in range(nchunks):
            lib_id = num + 1
            if lib_id not in library_info:
                library_info.update({lib_id: {}})
            val = aggr_df.iloc[num][key]
            if key in whitelist:
                if (key, val) in libraries:
                    return None, None, "Same value for {} is specified on multiple rows: {}".format(key, val)
                libraries.add((key, val))
            library_info[lib_id][key] = val
    return nchunks, library_info, None

def write_dict_to_csv(fn, dictionary, sort=False):
    """Basic dictionary to csv writer"""
    with open(fn, 'w') as f:
        keys = dictionary.keys() if not sort else sorted(dictionary.keys())
        f.write(','.join([str(k) for k in keys]) + '\n')
        f.write(','.join([str(dictionary[k]) for k in keys]) + '\n')


Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'strand'])
def peak_reader(peaks, select=[]):
    '''Streaming peak reader with selection mode. Accepts ordered indices'''
    curr = 0
    with open(peaks, 'r') as peak_file:
        for line, row in enumerate(peak_file):
            peak = row.strip("\n").split("\t")
            if curr < len(select) and line == select[curr]:
                curr += 1
                continue
            yield Peak(peak[0], int(peak[1]), int(peak[2]), '.')


def open_fragment_file(filename):
    with open_maybe_gzip(filename, "r") as fragment_file:
        for line in fragment_file:
            contig, start, stop, barcode, dup_count = line.strip().split("\t")
            yield contig, int(start), int(stop), barcode, int(dup_count)


def parsed_fragments_from_contig(contig, filename, index=None):
    """Returns an iterator yielding parsed fragments from a contig.
    Expects an index file of the same name as the input file, otherwise supply the file"""

    fragments = pysam.TabixFile(filename, index=index)
    try:
        for fragment in fragments.fetch(contig):
            contig, start, stop, barcode, dup_count = fragment.split("\t")
            yield contig, int(start), int(stop), barcode, int(dup_count)
    except ValueError:  # in the case when no regions on that contig are present
        pass

def grouped_fragments_from_contig(contig, filename, index=None):
    """Iterates through fragments, grouping up sets of fragments that share the same
    start position.
    """
    current_start = None
    fragment_list = []
    for fragment in parsed_fragments_from_contig(contig, filename, index):
        start = fragment[1]
        if current_start is None:
            current_start = start
        if start == current_start:
            fragment_list.append(fragment)
        else:
            yield fragment_list
            current_start = start
            fragment_list = [fragment]
    if fragment_list:
        yield fragment_list


def load_insert_sizes(filename, colrange=None, barcode_selection=None, report_sum=False):
    """Memory-efficient loading of the insert_size.csv file.  Uses chunking and loads in as uint32s.

    Optional arguments, all of which improve memory performance:
    - colrange: numerical values indicating columns (insert sizes) to return
    - barcode_selection: set of barcodes to report
    - report_sum: boolean if we want to report out the sum (excluding columns) instead of the dataframe.
    """
    with open(filename, "r") as infile:
        columns = infile.readline().strip().split(",")
    if colrange is None:
        usecols = None
    else:
        validcols = ["Barcode"] + [str(c) for c in colrange]

        def usecols(col):
            return col in validcols

    dtypes = {col: "string" if col == "Barcode" else "uint32" for col in columns}
    reader = pd.read_csv(filename, usecols=usecols, dtype=dtypes, chunksize=10000)
    chunks = []
    for chunk in reader:
        if barcode_selection is not None:
            cell_mask = chunk.Barcode.apply(lambda bc: bc in barcode_selection)
            chunk = chunk[cell_mask]
        if report_sum:
            # Exclude the barcode from the sum
            chunk = chunk.iloc[:, 1:].sum()
        chunks.append(chunk)

    if report_sum:
        return sum(chunks)
    else:
        return pd.concat(chunks)


class SubprocessStream(object):
    """Wrap a subprocess that we stream from or stream to.  Acts like an open filehandle by passing down
    next, fileno, write, and close down to its pipe.
    """

    def __init__(self, *args, **kwargs):
        mode = kwargs.pop("mode", "r")
        if mode == "r":
            kwargs["stdout"] = subprocess.PIPE
        elif mode == "w":
            kwargs["stdin"] = subprocess.PIPE
        else:
            raise ValueError("mode %s unsupported" % self.mode)

        kwargs["preexec_fn"] = os.setsid
        print args[0]
        sys.stdout.flush()
        self.proc = tk_subproc.Popen(*args, **kwargs)

        if mode == "r":
            self.pipe = self.proc.stdout
        elif mode == "w":
            self.pipe = self.proc.stdin

    def __enter__(self):
        return self

    def __iter__(self):
        return self

    def next(self):
        return self.pipe.next()

    def fileno(self):
        return self.pipe.fileno()

    def write(self, x):
        self.pipe.write(x)

    def close(self):
        self.pipe.close()
        self.proc.wait()

    def __exit__(self, tp, val, tb):
        self.close()


def open_maybe_gzip(filename, mode="r"):
    """Returns an open filehandle to a file that may or may not be compressed with gzip or lz4.
    """
    compressor = None

    if filename.endswith(GZIP_SUFFIX):
        compressor = "gzip"
    elif filename.endswith(LZ4_SUFFIX):
        compressor = "lz4"
    else:
        return open(filename, mode)

    assert compressor is not None

    if mode == "r":
        return SubprocessStream([compressor, "-c", "-d", filename], mode="r")
    elif mode == "w":
        f = open(filename, "w")
        return SubprocessStream([compressor, "-c"], stdout=f, mode="w")
    else:
        raise ValueError("Unsupported mode for compression: %s" % mode)


def combine_csv(input_csvs, output_csv, header_lines=1):
    """Combine a list of CSV files specified in input_csvs into
    a single output csv in output_csv. It is assumed that all
    CSVs have the same header structure. The number of header
    lines is specified in header_lines
    """
    if not input_csvs:
        output_csv = None
        return
    with open(output_csv, "w") as out:
        for i, icsv in enumerate(input_csvs):
            with open(icsv, "r") as infile:
                header = []
                for h in xrange(header_lines):
                    header.append(infile.next())
                if i == 0:
                    for hline in header:
                        out.write(hline)
                for line in infile:
                    out.write(line)


def create_bam_infile(file_name):
    # Note we set check_sq to False to allow correct handling of empty BAM files.
    bam_file = pysam.Samfile(file_name, "rb", check_sq=False)
    return bam_file


def hierarchical_merge_bam(input_bams, output_bam, tag=None, threads=1):
    """Merge sorted bam chunks hierarchically to conserve open file handles."""
    max_filehandles, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    max_filehandles = max(2, max_filehandles - 100)
    if len(input_bams) <= max_filehandles:
        # Special case this to avoid the shutil performance cost
        merge_bam(input_bams, output_bam, tag, threads)
    else:
        tmp_dir = os.path.dirname(output_bam)
        while len(input_bams) > 1:
            new_bams = []
            for i in range(0, len(input_bams), max_filehandles):
                bam_chunk = input_bams[i:i + max_filehandles]
                if len(bam_chunk) > 1:
                    new_bam = os.path.join(tmp_dir, "%d-%d.bam" % (i, len(input_bams)))
                    merge_bam(bam_chunk, new_bam, tag, threads)
                    new_bams.append(new_bam)
                else:
                    new_bams.append(input_bams[i])
            input_bams = new_bams
        shutil.move(input_bams[0], output_bam)


def merge_bam(input_bams, output_bam, tag=None, threads=1):
    """ Merge all input files into a single output file, keeping sorted order on a given tag
    """
    args = []
    if tag is not None:
        args.extend(["-t", tag])
    if threads > 1:
        args.extend(["-c", "-p", "-s", "0", "-@", threads])
    args.append(output_bam)
    args.extend(list(input_bams))
    args = [str(arg) for arg in args]
    pysam.merge(*args)


def sort_bam(input_bam, output_bam, tag=None, threads=1):
    args = ["-o", output_bam]
    if tag is not None:
        args.extend(["-t", tag])
    if threads > 1:
        args.extend(["-@", threads])
    args.append(input_bam)
    args = [str(arg) for arg in args]
    pysam.sort(*args)


def index_bam(input_bam, threads=1):
    args = [input_bam]
    if threads > 1:
        args.extend(["-@", threads])
    args = [str(arg) for arg in args]
    pysam.index(*args)
    index_file = "{}.bai".format(input_bam)
    if not os.path.isfile(index_file):
        raise RuntimeError("samtools index failed, likely bam is not sorted")
    return index_file


def sort_bed(input_bed, output_bed, genome, threads=1, leave_key=False, has_key=False):
    """Use unix sort to properly sort a bed file, including a custom sort order on chromosomes.

    Warning!  sort does not have the --parallel argument on all forms of unix!  As such, we are
    dropping threading support for BED handling through unix sort.
    """
    main_cmds = []
    # sort_thread_args = "" if threads == 1 else " --parallel={}".format(threads)
    sort_thread_args = ""
    if not has_key:
        # If the bed file doesn't already have a contig key, we need to create the key file used to add it.
        tmpdir = os.path.dirname(output_bed)
        tmp_chroms = os.path.join(tmpdir, "chrom_order.txt")
        chroms = os.path.join(tmpdir, "sorted_order.txt")
        with open(genome, "r") as infile, open(tmp_chroms, "w") as outfile:
            for i, line in enumerate(infile):
                chrom = line.split()[0]
                outfile.write("{}\t{}\n".format(chrom, i))
        tk_subproc.check_call("sort -k1b,1 -o {} {}".format(chroms, tmp_chroms), shell=True)
        # Now we'll add the commands to join the key onto the file between the contig & start/stop positions
        main_cmds.extend(["sort -k1b,1{} {}".format(sort_thread_args, input_bed),
                          "join -t '\t' -j1 {} -".format(chroms)])
    else:
        main_cmds.append("cat {}".format(input_bed))
    # Next we sort on the contig key and start/stop positions
    main_cmds.append("sort -k2n -k3n -k4n{}".format(sort_thread_args))

    if not leave_key:
        # Finally we remove the key from the output file
        main_cmds.append("cut -f 1,3-8")

    with open(output_bed, 'w') as outfile:
        tk_subproc.check_call(" | ".join(main_cmds), shell=True, stdout=outfile)


def merge_keyed_bed(input_beds, output_bed, threads=1):
    """Merge sorted bedfiles retaining their chromosome keys, dropping the key afterwards.

    Warning!  sort does not have the --parallel argument on all forms of unix!  As such, we are
    dropping threading support for BED handling through unix sort.
    """
    main_cmds = []
    # sort_thread_args = "" if threads == 1 else " --parallel={}".format(threads)
    sort_thread_args = ""
    main_cmds.append("sort -m -k2n -k3n -k4n{} {}".format(sort_thread_args, ' '.join(input_beds)))
    main_cmds.append("cut -f 1,3-8")

    with open(output_bed, 'w') as outfile:
        tk_subproc.check_call(" | ".join(main_cmds), shell=True, stdout=outfile)

def subsample_fragments(infile, rate, outfile, cells=set(), group=None, kind=None, outindex=None, key=None, randomseed=0):
    """Downsample a fragments infile to outfile (and optionally build the tabindex). Downsampling can be of 3 kind:
    - None: simply copies the input
    - unique: downsamples the number of unique fragments
    - total: downsample the total number of fragments
    If provided cell calls, it computes single cell metrics pre and post downsampling
    If provided group, it relabels the barcodes to be from that gem group
    If provided outindex, it computes a tabix and compresses the downsampled file in process
    If provided key, it adds a extra column about chromosome order, to be used in merging multiple files based on a specific order
    Different bulk library sensitivity metrics pre and post normalization are calculated
    """
    assert kind in [None, "total", "unique"]
    metrics = {'total_pre_normalization': 0,
               'total_post_normalization': 0,
               'unique_pre_normalization': 0,
               'unique_post_normalization': 0}

    np.random.seed(randomseed)

    basefile, _ = os.path.splitext(outfile)
    final_outfile = basefile if outindex is not None else outfile

    if group is not None:
        gem_group_suff = "-{}".format(group)

    np.random.seed(0)

    pre_frags_per_cell = Counter({bc: 0 for bc in cells})
    post_frags_per_cell = Counter()
    with open(final_outfile, 'w') as f:
        for contig, start, stop, barcode, dup_count in open_fragment_file(infile):
            is_cell = False
            if barcode in cells:
                is_cell = True
            if is_cell:
                pre_frags_per_cell[barcode] += 1
            if group is not None:
                barcode = barcode.strip("-1") + gem_group_suff
            metrics['total_pre_normalization'] += dup_count
            metrics['unique_pre_normalization'] += 1
            if kind == "total":
                new_dup_count = sum(np.random.uniform(size=dup_count) < rate)
                if new_dup_count > 0:
                    if key is None:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(contig, start, stop, barcode, new_dup_count))
                    else:
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, key[contig], start, stop, barcode, new_dup_count))
                    metrics['total_post_normalization'] += new_dup_count
                    metrics['unique_post_normalization'] += 1
                    if is_cell:
                        post_frags_per_cell[barcode] += 1
            if kind == "unique" or kind is None:
                if kind is None or np.random.uniform(size=1) < rate:
                    if key is None:
                        f.write("{}\t{}\t{}\t{}\t{}\n".format(contig, start, stop, barcode, dup_count))
                    else:
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, key[contig], start, stop, barcode, dup_count))
                    metrics['total_post_normalization'] += dup_count
                    metrics['unique_post_normalization'] += 1
                    if is_cell:
                        post_frags_per_cell[barcode] += 1
    if len(cells) > 0:
        metrics['pre_norm_median_frags_per_cell'] = np.median(pre_frags_per_cell.values())
        metrics['post_norm_median_frags_per_cell'] = np.median(post_frags_per_cell.values())

    # convert to gz and tabix
    if outindex is not None:
        pysam.tabix_index(basefile, preset='bed', index=outindex)
    return metrics


def get_counts_by_barcode(reference_path, peaks, fragments, fragments_index=None, contig=None, known_cells=None):
    """Generate targeting, raw and dup counts per barcode. If cell identity is known, then also return that as part of
    the counts
    """
    def load_reference_track(track, padding=0):
        if track is not None:
            with open(track, 'r') as infile:
                regions = regtools.get_target_regions(infile, padding=padding)
        else:
            regions = None
        return regions

    def point_is_in_target(contig, position, target_regions):
        if target_regions is None:
            return False
        if contig not in target_regions:
            return False
        return target_regions[contig].contains_point(position)

    def fragment_overlaps_target(contig, start, stop, target_regions):
        if target_regions is None:
            return False
        if contig not in target_regions:
            return False
        return target_regions[contig].overlaps_region(start, stop)

    ref_manager = ReferenceManager(reference_path)

    # Load in and pad TSS/CTCF regions if present
    tss_regions = load_reference_track(ref_manager.tss_track, padding=2000)
    ctcf_regions = load_reference_track(ref_manager.ctcf_track, padding=250)

    # Load in regions from reference-associated tracks
    dnase_regions = load_reference_track(ref_manager.dnase_track)
    enhancer_regions = load_reference_track(ref_manager.enhancer_track)
    promoter_regions = load_reference_track(ref_manager.promoter_track)
    blacklist_regions = load_reference_track(ref_manager.blacklist_track)
    peak_regions = load_reference_track(peaks)

    # load cell - species map
    cell_barcodes = {}
    species_list = ref_manager.list_species()
    if known_cells is not None:
        with open(known_cells, 'r') as infile:
            for line in infile:
                items = line.strip("\n").split(",")
                for barcode in items[1:]:
                    if barcode != "null":
                        if barcode not in cell_barcodes:
                            cell_barcodes[barcode] = []
                        cell_barcodes[barcode] += [items[0]]

    # get cell index
    cell_index = {}
    spnum = {species: 0 for species in species_list}
    for species in species_list:
        for barcode in cell_barcodes:
            if species in cell_barcodes[barcode]:
                label = "{}_cell_{}".format(species, spnum[species])
                spnum[species] += 1
                cell_index[barcode] = label if barcode not in cell_index else '_'.join([cell_index[barcode], label])

    counts_by_barcode = {}
    tss_relpos = Counter()
    ctcf_relpos = Counter()

    read_count = 0

    iterator = open_fragment_file(fragments) if contig is None else \
        parsed_fragments_from_contig(contig, fragments, index=fragments_index)
    for contig, start, stop, barcode, dups in iterator:
        read_count += 2
        if barcode not in counts_by_barcode:
            counts_by_barcode[barcode] = Counter()
            if known_cells is not None:
                cell_species = cell_barcodes.get(barcode, [])
                counts_by_barcode[barcode]["cell_id"] = cell_index.get(barcode, "None")
                for species in species_list:
                    if species in cell_species:
                        counts_by_barcode[barcode]["is_{}_cell_barcode".format(species)] = 1
                    else:
                        counts_by_barcode[barcode]["is_{}_cell_barcode".format(species)] = 0

        # species splits
        if known_cells is not None and len(species_list) > 1:
            contig_species = ref_manager.species_from_contig(contig)
            counts_by_barcode[barcode]["passed_filters_{}".format(contig_species)] += 1
            if fragment_overlaps_target(contig, start, stop, peak_regions):
                counts_by_barcode[barcode]["peak_region_fragments_{}".format(contig_species)] += 1

        # raw mapping
        counts_by_barcode[barcode]["passed_filters"] += 1
        counts_by_barcode[barcode]["total"] += dups
        counts_by_barcode[barcode]["duplicate"] += dups - 1

        # Count up transposition site targeting
        for position in (start, stop):
            if point_is_in_target(contig, position, tss_regions):
                region = tss_regions[contig].get_region_containing_point(position)
                tss_relpos[region.get_relative_position(position)] += 1
            if point_is_in_target(contig, position, ctcf_regions):
                region = ctcf_regions[contig].get_region_containing_point(position)
                ctcf_relpos[region.get_relative_position(position)] += 1
            if point_is_in_target(contig, position, peak_regions):
                counts_by_barcode[barcode]["peak_region_cutsites"] += 1

        # Count up fragment overlap targeting
        is_targeted = False
        if fragment_overlaps_target(contig, start, stop, tss_regions):
            counts_by_barcode[barcode]["TSS_fragments"] += 1
            is_targeted = True
        if fragment_overlaps_target(contig, start, stop, dnase_regions):
            counts_by_barcode[barcode]["DNase_sensitive_region_fragments"] += 1
            is_targeted = True
        if fragment_overlaps_target(contig, start, stop, enhancer_regions):
            counts_by_barcode[barcode]["enhancer_region_fragments"] += 1
            is_targeted = True
        if fragment_overlaps_target(contig, start, stop, promoter_regions):
            counts_by_barcode[barcode]["promoter_region_fragments"] += 1
            is_targeted = True
        if is_targeted:
            counts_by_barcode[barcode]["on_target_fragments"] += 1
        if fragment_overlaps_target(contig, start, stop, blacklist_regions):
            counts_by_barcode[barcode]["blacklist_region_fragments"] += 1
        if fragment_overlaps_target(contig, start, stop, peak_regions):
            counts_by_barcode[barcode]["peak_region_fragments"] += 1
    return read_count, counts_by_barcode, tss_relpos, ctcf_relpos
