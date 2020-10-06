#!/usr/bin/env python
#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#

from __future__ import division

import sys
import urllib
import re
import os
import json
import subprocess
import shutil
import requests
import numpy as np
import pyfasta

from tools import ReferenceManager
from tools.io import peak_reader
from sklearn.utils import sparsefuncs
import pandas as pd

def normalize_matrix(matr, scale):
    """normalize a matrix with some scale"""
    m = matr.copy().astype(np.float64)
    scale = np.median(scale) / scale
    sparsefuncs.inplace_column_scale(m, scale)
    return m

def get_barcode_gc(ref_f, peaks_f, matrix):
    """Get mean GC% of peaks in a barcode"""
    ref_mgr = ReferenceManager(ref_f)
    genome_fa = pyfasta.Fasta(ref_mgr.fasta, key_fn=lambda x: x.split()[0])
    peak_GC = np.array([get_peak_GC_counts(peak, genome_fa, counts=False)
                        for peak in peak_reader(peaks_f)])
    barcode_GC = ((peak_GC * matrix.m) / np.array(matrix.m.sum(axis=0))).squeeze()
    return barcode_GC

def load_cell_barcodes(filename, with_species=False):
    """Read cell_barcodes.csv and emit barcodes"""
    cell_barcodes = {}
    with open(filename, 'r') as f:
        for line in f:
            items = line.strip("\n").rstrip(",").split(",")
            cell_barcodes[items[0]] = set({item for item in items[1:] if item != "null" and item != ""}) if len(items) > 1 else set()
    if not with_species:
        bcs = set()
        for species in cell_barcodes:
            bcs.update(cell_barcodes[species])
        return bcs
    else:
        return cell_barcodes

def get_cell_barcodes(filename, ref, with_species=False):
    """Read singlecell.csv and emit barcodes"""
    scdf = pd.read_csv(filename, sep=',')
    ctg_mgr = ReferenceManager(ref)
    if not with_species:
        cell_barcodes = set()
        for species in ctg_mgr.list_species():
            species_cell_mask = scdf['is_{}_cell_barcode'.format(species)] == 1
            cell_barcodes.update(scdf[species_cell_mask]['barcode'].values.tolist())
    else:
        cell_barcodes = {}
        for species in ctg_mgr.list_species():
            species_cell_mask = scdf['is_{}_cell_barcode'.format(species)] == 1
            cell_barcodes[species] = set(scdf[species_cell_mask]['barcode'].values.tolist())
    return cell_barcodes

def get_peak_GC_counts(peak, genome_fa, counts=True):
    '''Get GC% in base seq in a peak and (optionally) nucleotide counts'''

    seq = genome_fa[peak.chrom][peak.start:peak.end].upper()
    base_counter = {}
    base_counter['A'] = seq.count('A') + seq.count('a')
    base_counter['G'] = seq.count('G') + seq.count('g')
    base_counter['C'] = seq.count('C') + seq.count('c')
    base_counter['T'] = seq.count('T') + seq.count('t')
    base_counter['N'] = seq.count('N') + seq.count('n')
    peakGC = (base_counter['G'] + base_counter['C']) / sum(base_counter.values()) if len(seq) > 0 else 0.
    if counts:
        return peakGC, base_counter
    else:
        return peakGC

def generate_genome_tag(ref_path):
    """Replace empty genome name for single genomes with valid genome name"""
    # For a single species reference, use contents of <reference_path>/genome
    ref_contig_manager = ReferenceManager(ref_path)
    genomes = ref_contig_manager.list_species()
    if len(genomes) == 1 and genomes[0] == '' or len(genomes) == 0:
        genomes = [ref_contig_manager.genome]
    return genomes

def count_bases_in_peaks(reference_path, peaks_file):
    """Count the total number of bases in peak regions (0-indexed)"""
    bases_in_peaks = 0
    ctg_mgr = ReferenceManager(reference_path)
    genome_fa = pyfasta.Fasta(ctg_mgr.fasta, key_fn=lambda x: x.split()[0])
    for peak in peak_reader(peaks_file):
        bases_in_peaks += len(genome_fa[peak.chrom][peak.start:peak.end])
    return bases_in_peaks

def cut_site_counter_from_bedgraph(cutsites_f, threshold=1, contig=None):
    """Yield count data at base pairs encoded in bedgraph tracks"""
    with open(cutsites_f, 'r') as f:
        for track in f:
            chrom, start, stop, count = track.split("\t")
            if contig is not None and chrom != contig:
                continue
            if count < threshold:
                continue
            for pos in xrange(int(start), int(stop)):
                yield chrom, pos, int(count)

def get_purity_info(singlecell_df, species_list):
    '''Calculate species barcode purity'''
    assert len(species_list) == 2

    def calculate_purity(self_counts, other_counts):
        return self_counts / (self_counts + other_counts)

    xdata = singlecell_df["passed_filters_{}".format(species_list[0])].values
    ydata = singlecell_df["passed_filters_{}".format(species_list[1])].values
    is_spec_1 = singlecell_df['is_{}_cell_barcode'.format(species_list[0])] == 1
    is_spec_2 = singlecell_df['is_{}_cell_barcode'.format(species_list[1])] == 1
    spec_1_purity = calculate_purity(xdata[is_spec_1], ydata[is_spec_1])
    spec_2_purity = calculate_purity(ydata[is_spec_2], xdata[is_spec_2])

    return spec_1_purity, spec_2_purity

def get_chunks(num, chunks):
    '''Attempts to generate non-empty slices out of num entries.'''
    bounds = np.linspace(0, num, chunks + 1, dtype=int, endpoint=True).tolist()
    return [(s, e) for s, e in zip(bounds, bounds[1:]) if s != e]

def download_file(url, target_path):
    if url.startswith("ftp"):
        urllib.urlretrieve(url, filename=target_path)
        return

    r = requests.get(url, stream=True)
    with open(target_path, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

def download_from_url(url, target_path):

    if url.endswith('.gz'):
        target_path += '.gz'
        download_file(url, target_path)
        subprocess.call(['gzip', '-f', '-d', target_path])

    else:
        download_file(url, target_path)

def is_url(input_url):
    # detect the input string is a valid url, borrowed from django URLValidator

    ul = '\u00a1-\uffff'  # unicode letters range (must not be a raw string)

    # IP patterns
    ipv4_re = r'(?:25[0-5]|2[0-4]\d|[0-1]?\d?\d)(?:\.(?:25[0-5]|2[0-4]\d|[0-1]?\d?\d)){3}'
    ipv6_re = r'\[[0-9a-f:\.]+\]'  # (simple regex, validated later)

    # Host patterns
    hostname_re = r'[a-z' + ul + r'0-9](?:[a-z' + ul + r'0-9-]{0,61}[a-z' + ul + r'0-9])?'
    # Max length for domain name labels is 63 characters per RFC 1034 sec. 3.1
    domain_re = r'(?:\.(?!-)[a-z' + ul + r'0-9-]{1,63}(?<!-))*'
    tld_re = (
        r'\.'  # dot
        r'(?!-)'  # can't start with a dash
        r'(?:[a-z' + ul + '-]{2,63}'  # domain label
        r'|xn--[a-z0-9]{1,59})'  # or punycode label
        r'(?<!-)'  # can't end with a dash
        r'\.?'  # may have a trailing dot
    )

    # Note: modified from original, allows tld (top level domain name) missing
    # host_re = '(' + hostname_re + domain_re + tld_re + '|localhost)'
    host_re = '(' + hostname_re + domain_re + '|' + tld_re + '|localhost)'

    valid_url_regex = re.compile(
        r'^(?:[a-z0-9\.\-\+]*)://'  # scheme is validated separately
        r'(?:[^\s:@/]+(?::[^\s:@/]*)?@)?'  # user:pass authentication
        r'(?:' + ipv4_re + '|' + ipv6_re + '|' + host_re + ')'
        r'(?::\d{2,5})?'  # port
        r'(?:[/?#][^\s]*)?'  # resource path
        r'\Z', re.IGNORECASE)

    return re.match(valid_url_regex, input_url) is not None

def validate_file(input_file, filetype):
    # If an input file is provided, verify that the input file exists

    if not os.path.exists(input_file):
        sys.exit("Input {} file does not exist: {}".format(filetype, input_file))

    if not os.path.isfile(input_file):
        sys.exit("Please provide a file, not a directory: {}".format(input_file))

def fetch_from_source(source, target_path, filetype):

    # If the targeted file already exists, skip the copy/download step
    if os.path.exists(target_path):
        print("{} file already exists! Skip fetching from source...\n".format(target_path))

    else:
        # If the input source is an url, download it. Otherwise treat it as a path
        if is_url(source):
            print "Downloading {} files from source...".format(filetype)
            download_from_url(source, target_path)

        else:
            # Copy source file to target file
            print "Copying original {} file into reference folder...".format(filetype),

            # valid input file, if failed, program will exit with an error message
            validate_file(source, filetype)
            if source.endswith('.gz'):
                subprocess.call(['gzip', '-c', '-d', source], stdout=open(target_path, 'w'))

            else:
                shutil.copy(source, target_path)
            os.chmod(target_path, 0644)
        print "done\n"

def extact_contig_names_from_def(ref_path):
    contig_defs_ready_file = os.path.join(ref_path, 'fasta', 'contig-defs.json')

    contig_defs_ready = json.load(open(contig_defs_ready_file))
    contig_names = contig_defs_ready['primary_contigs'] + contig_defs_ready['non_nuclear_contigs']

    return contig_names


def quick_line_count(filename):
    """
    a fast line counter. It reads the text file as binary chunks and count the number of new line characters
    :param filename: str of path of the input file
    :return: number of lines
    """
    def _make_gen(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024 * 1024)

    f = open(filename, 'rb')
    f_gen = _make_gen(f.read)
    return sum(buf.count(b'\n') for buf in f_gen)

# check format in motifs.pfm
def motif_format_checker(motif_infile):
    from Bio import motifs

    if motif_infile is None:
        return None

    try:
        all_motifs = motifs.parse(open(motif_infile), "jaspar")
    except:
        sys.exit("Motif file is not in JASPAR format.")

    nmotif = 0
    with open(motif_infile) as motif_in:
        for i, line in enumerate(motif_in):
            if line.startswith('>'):
                nmotif += 1
                if len(line.split('\t')) > 1:
                    sys.exit(
                        "Motif name cannot contain tabs('\t') at line {} in {}.".format(i+1, motif_infile)
                    )

    if nmotif != len(all_motifs):
        sys.exit("Motif file is not in JASPAR format.")
