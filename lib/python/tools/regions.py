"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Data structure and functions for region-based file manipulation, such ad BED, GFF and GTF.

Originally copied from tenkit and simplified.
"""

import bisect
import json
import os
import shutil
import subprocess
import sys
import tempfile

from collections import namedtuple
from pybedtools import BedTool, MalformedBedLineError


class Region(namedtuple('Region', ['start', 'end', 'name', 'strand'])):
    """Overloads a namedtuple to add a few additional methods.
    """
    @property
    def width(self):
        return self.end - self.start - 1

    @property
    def is_positive_strand(self):
        return not self.strand == '-'

    def get_relative_position(self, point):
        if self.is_positive_strand:
            return (point - self.start) - self.width // 2
        else:
            return (self.end - point) - self.width // 2


class Regions(object):
    """Defines a set of regions and allows for determining whether a point
    is contained in any of those regions or whether another region overlaps any of them.
    """

    def __init__(self, region_list=None):
        if region_list is not None and not (region_list == []):
            sorted_regions = sorted(region_list)
            self.starts = [r[0] for r in sorted_regions]
            self.ends = [r[1] for r in sorted_regions]
            self.region_list = sorted_regions
        else:
            self.starts = []
            self.ends = []
            self.region_list = []
        self.current_iter = 0

    def __iter__(self):
        return self

    def next(self):
        if self.current_iter >= len(self.region_list):
            self.current_iter = 0
            raise StopIteration
        else:
            self.current_iter += 1
            return self.region_list[self.current_iter]

    def contains_point(self, point):
        """ Determines whether a point is contained in one of the regions.
        """
        index = bisect.bisect(self.starts, point) - 1
        return index >= 0 and point >= self.starts[index] and point < self.ends[index]

    def get_region_containing_point(self, point):
        """ Returns a region containing the point, if one exists.
        """
        index = bisect.bisect(self.starts, point) - 1
        if index == -1:
            return None

        if point >= self.starts[index] and point < self.ends[index]:
            return self.region_list[index]
        return None

    def overlaps_region(self, start, end):
        """ Determines whether a region described by start and end overlaps
        any of the regions.
        """
        check_index = bisect.bisect_left(self.ends, start)
        if check_index == len(self.starts):
            return False
        if end > self.starts[check_index]:
            return True
        return False


def get_target_regions_dict(targets_file, padding=0):
    """ Gets the target regions from a targets BED file as a chrom-indexed dictionary,
    with every entry given as a list of (start, end, name, strand) tuples
    """
    targets = {}
    for line in targets_file:
        info = line.strip().split('\t')
        if any(line.startswith(comment) for comment in ['browser', 'track', '-browser', '-track', '#']):
            continue
        if len(line.strip()) == 0:
            continue

        chrom = info[0]
        start = int(info[1]) - padding
        end = int(info[2]) + padding
        name = info[3] if len(info) >= 4 else None
        strand = info[5] if len(info) >= 6 else None
        if strand not in ['+', '-']:
            strand = None
        chrom_targs = targets.setdefault(chrom, [])
        chrom_targs.append(Region(start, end, name, strand))
    return targets


def get_target_regions(targets_file, padding=0):
    """ Gets the target regions from a targets file as a chrom-indexed dictionary,
    with every entry given as a list of Regions objects
    """
    targets_dict = get_target_regions_dict(targets_file, padding)
    target_regions = {}
    for (chrom, region_list) in targets_dict.iteritems():
        target_regions[chrom] = Regions(region_list=region_list)
    return target_regions


def fragment_overlaps_target(contig, start, stop, target_regions):
    if target_regions is None:
        return False
    if contig not in target_regions:
        return False
    return target_regions[contig].overlaps_region(start, stop)


def clean_chr_name(chrom, contig_defs, keep_nonstd=False):
    """
    match the input chrom string to the standard format used in the genome.fa
    """

    std_names = contig_defs['primary_contigs'] + contig_defs['non_nuclear_contigs']
    test_name = chrom.lstrip('chr')
    if test_name in std_names:
        return contig_defs['prefix'] + test_name

    # edge cases for mitochondrial chrom
    elif test_name == 'M':
        return contig_defs['prefix'] + 'MT'
    elif test_name == 'MT':
        return contig_defs['prefix'] + 'M'
    elif keep_nonstd:
        return test_name
    else:
        return None


def clean_chr_name_file(in_file, contig_defs, ref_path):
    """
    given a tab-delimitated file, run clean_chr_name to clean the first column, which is the default chrom column
    :return name changes will be make in the original file
    """

    print("Filtering {}...".format(in_file))
    # run validate_contig_names first before running chr name clean. if pass, no clean_chr_name is needed
    if ref_path is not None:
        contig_defs_ready_file = os.path.join(ref_path, 'fasta', 'contig-defs.json')
        if os.path.exists(contig_defs_ready_file):
            with open(contig_defs_ready_file, 'r') as contig_defs_ready_file_obj:
                contig_defs_ready = json.load(contig_defs_ready_file_obj)
            contig_defs_ready['PRIMARY_CONTIGS'] = contig_defs_ready['primary_contigs']
            contig_defs_ready['NON_NUCLEAR_CONTIGS'] = contig_defs_ready['non_nuclear_contigs']
            try:
                validate_contig_names(in_file, contig_defs_ready)
                return
            except:
                pass

    tmp = tempfile.NamedTemporaryFile(dir=os.path.dirname(in_file), delete=False)

    bediter = bed_format_chcker_iter(in_file)

    for fields in bediter:
        chrom = clean_chr_name(str(fields.chrom), contig_defs)

        if chrom is not None:
            fields.chrom = chrom
            tmp.write('\t'.join(fields.fields) + '\n')

    # move the processed file to the input file
    shutil.move(tmp.name, in_file)
    tmp.close()
    print("done\n")


def gtf_or_gff3(in_file, check_first_n=100):
    """
    Check an input file is a GTF file or a GFF3 file
    :param in_file: input file name
    :param check_first_n: check the first n rows of the file. Default = 100
    :return:
    """
    gxf_type = {'GTF': 0, 'GFF': 0}
    bediter = bed_format_chcker_iter(in_file)
    i = count_header_lines_bed(in_file)
    check_first_n += i
    for row in bediter:
        i += 1
        first_attr_str = row.fields[8].split(';')[0]
        first_attr_key = first_attr_str.split(' ')[0]

        try:
            row.attrs
        except ValueError:
            msg = "Malformed GTF/GFF attributes at line {} in {}".format(i, in_file)
            sys.exit(msg)

        if first_attr_key in row.attrs:
            gxf_type['GTF'] += 1

        first_attr_key = first_attr_str.split('=')[0]
        if first_attr_key in row.attrs:
            gxf_type['GFF'] += 1

        if i > check_first_n:
            break

    if gxf_type['GTF'] > gxf_type['GFF'] and gxf_type['GFF'] == 0:
        return 'GTF'

    elif gxf_type['GTF'] < gxf_type['GFF'] and gxf_type['GTF'] == 0:
        return 'GFF'

    else:
        sys.exit("Input file is neither a valid GTF nor a valid GFF file: {}".format(in_file))


def convert_gff3_to_gtf(gff_infile, gtf_outfile):
    """
    Convert GFF3 format file into GTF format file
    :param gff_infile: input file path
    :param gtf_outfile: output file path
    """

    def _name_to_genename(attr_dict):
        """replace Name of GFF3 attributes key with gene_name"""

        new_attr_dict = dict(attr_dict)

        if "gene_name" not in attr_dict:
            if "Name" in attr_dict:
                new_attr_dict["gene_name"] = attr_dict["Name"]

            # parse Parent when gene_name and Name are both absent
            # covers format as Parent=gene:XXX or Parent=XXX by taking the last element
            # of the value splitted by colon
            elif "Parent" in attr_dict:
                new_attr_dict["gene_name"] = attr_dict["Parent"].split(":")[-1]

        return new_attr_dict

    new_gtf = []
    new_gtf_writer = open(gtf_outfile, 'w')

    gtf_iter = bed_format_chcker_iter(gff_infile)
    i = count_header_lines_bed(gff_infile)
    for row in gtf_iter:
        i += 1

        try:
            new_attrs = _name_to_genename(row.attrs)
        except ValueError:
            msg = "Malformed GTF/GFF attributes at line {} in {}".format(i, gff_infile)
            sys.exit(msg)

        new_attr_fields = [
            '{} "{}"'.format(k, v) for k, v in new_attrs.iteritems()
        ]
        new_attr_field = ";".join(new_attr_fields) + ";"
        new_row = "\t".join(row.fields[:8] + [new_attr_field]) + "\n"
        new_gtf.append(new_row)

        # write 100k lines at a time
        if i % 100000 == 0:
            new_gtf_writer.writelines(new_gtf)
            new_gtf = []

    new_gtf_writer.writelines(new_gtf)
    new_gtf_writer.close()

    assert gtf_or_gff3(gtf_outfile) == 'GTF'


def is_overlapping(in_peaks_file):
    """Given a sorted bed file of peaks, check if the peaks are non overlapping"""

    assert in_peaks_file is not None
    last_chrom = None
    last_peak = None

    bediter = bed_format_chcker_iter(in_peaks_file)
    while bediter:
        try:
            bed = bediter.next()
            if bed.chrom != last_chrom:
                last_chrom = bed.chrom
                last_peak = (bed.start, bed.end)
            else:
                if bed.start < last_peak[1]:
                    return True
        except StopIteration:
            break
    return False


def bed_format_chcker_iter(in_bed_file):
    """
    Take a file path of a BED/GTF/GFF file and return a generator with BedTool built-in format checker
    :param in_bed_file:
    :return:
    """

    bediter = iter(BedTool(in_bed_file))
    i = count_header_lines_bed(in_bed_file)

    while True:
        i += 1
        try:
            row = bediter.next()
            yield row
        except MalformedBedLineError as e:
            e_is_too_long = '\n' if len(str(e)) > 80 else ''
            msg = "Malformed BED entry at line {} in {}: {}".format(i, in_bed_file, e_is_too_long)
            sys.exit(msg + str(e))
        except IndexError as e:
            if i == 1:
                msg = 'Invalid BED format for {}, please check the source.'.format(in_bed_file)
                sys.exit(msg)
            else:
                msg = "Malformed BED entry at line {} in {}: ".format(i, in_bed_file)
                sys.exit(msg + str(e))
        except StopIteration:
            raise


def bed_format_checker(in_bed_file, faidx_file):
    """
    Given a bed file path, examine whether is it a valid bed format
    Examined items include: basic formatting, coordinate overflow, contig name is within reference and properly sorted
    :param in_bed_file: file path of input bed file
    :param faidx_file: fasta index file (.fai) file.
    :return:
    """
    if in_bed_file is None:
        return None

    all_contigs = []
    contig_lengths = {}
    with open(faidx_file) as faidx:
        for line in faidx:
            line = line.strip().split()
            contig_name, contig_length = line[0], line[1]

            all_contigs.append(contig_name)
            contig_lengths[contig_name] = int(contig_length)

    bediter = bed_format_chcker_iter(in_bed_file)
    contig_iter = iter(all_contigs)  # contig order in the reference
    next_contig = next(contig_iter)

    i = count_header_lines_bed(in_bed_file)
    for bed in bediter:
        i += 1
        if bed.chrom not in contig_lengths:
            sys.exit('Malformed BED entry at line {} in {}: invalid contig name.'.format(i, in_bed_file))

        if bed.end > contig_lengths[bed.chrom]:
            if os.path.basename(in_bed_file) == 'blacklist.bed':
                pass
            else:
                sys.exit('Malformed BED entry at line {} in {}: Coordinate exceeds the length of the contig.'.format(i, in_bed_file))

        while bed.chrom != next_contig:
            try:
                next_contig = next(contig_iter)
            except StopIteration:
                sys.exit('Malformed BED entry at line {} in {}: contigs are not properly sorted according to the order in the reference fasta.'.format(i, in_bed_file))


def sort_and_uniq_bed(in_bed_file, ref_path):
    """
    given a path of a bed file, sort bed file and remove duplicate rows
    ref_path is required to extract the chromosome orders from the genome.fa.fai
    """

    if in_bed_file.endswith('.gz'):
        sys.exit("Input bed file is gzip compressed, please decompress it first.\n")

    # chrom are sorted numerically
    faidx_file = os.path.join(ref_path, 'fasta', 'genome.fa.fai')
    if not os.path.exists(faidx_file):
        sys.exit("Cannot find genome fasta index file.")

    input_bed = BedTool(in_bed_file)

    with tempfile.NamedTemporaryFile(dir=os.path.dirname(in_bed_file), delete=False) as tmp:
        sorted_bed = input_bed.sort(g=faidx_file).saveas(tmp.name)

        # remove duplicate lines and write to the original input bed file
        uniq_cmd = ['uniq', sorted_bed.fn]
        with open(in_bed_file, 'w') as out_file:
            subprocess.call(uniq_cmd, stdout=out_file)

        # clean up temporary files
        os.remove(tmp.name)

def validate_contig_names(in_file, contig_defs):
    """
    Validate contig names for fasta index or BED/GTF/GFF files.
    All contigs specified in "PRIMARY_CONTIGS" and "NON_NUCLEAR_CONTIGS" of <contig_defs> have to be present in <in_file>
    :param in_file: samtools faidx file
    :param contig_defs: a dict containing definition of PRIMARY_CONTIGS
    :return:
    """

    filetype = os.path.splitext(in_file)[1].lstrip('.')

    contigs_seen = set()
    contigs_in_def = set(contig_defs['PRIMARY_CONTIGS'] + contig_defs['NON_NUCLEAR_CONTIGS'])

    if filetype == 'fai':
        with open(in_file) as faidx:
            for line in faidx:
                contig_name = line.strip().split('\t')[0]
                contigs_seen.add(contig_name)

    elif filetype in ['bed', 'gtf', 'gff', 'gff3']:
        # extract the first column by cut and run uniq, faster than using BedTool and iterate each row
        contig_col = subprocess.Popen(["cut", "-f", "1", in_file], stdout=subprocess.PIPE)
        contig_col_uniq = subprocess.check_output(('uniq'), stdin=contig_col.stdout)
        contigs_seen = set(contig_col_uniq.strip().split('\n'))

    else:
        sys.exit('Invalid file extension for {}'.format(in_file))

    if len(contigs_in_def - contigs_seen) > 0:
        sys.exit("Invalid contig names. All PRIMARY_CONTIGS in the configuration file need to be present in {}". format(in_file))
    else:
        return

def count_header_lines_bed(in_file):
    """ Count the number of header lines in a bed file
        Header lines are defined as lines starting with "#" in the beginning of the file
        Comment lines with "#" that not in the header will not be counted
    """

    counter = 0
    with open(in_file) as file_in:
        for line in file_in:
            if line.startswith("#"):
                counter += 1
            else:
                break

    return counter