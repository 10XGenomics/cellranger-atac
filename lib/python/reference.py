"""
Source files for genome annotation. Each variable is a dictionary in which
keys are supported genome versions
    GTF_SOURCE: GENCODE basic CHR gtf file
    DHS_SOURCE: Master Peak List ENCODE3 human and mouse DNase-seq October 2015
    BLACKLIST_SOURCE:   ENCODE black list for ChIP-seq, DNase-sea and ATAC-seq
    REG_ENSEMBL:    Ensembl regulatory build

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

# TODO: Write up proper docstrings for this module
# pylint: disable=missing-docstring

from __future__ import division, print_function, absolute_import

import collections
import csv
import json
import os
import shutil
import subprocess
import sys

import hjson
import pyfasta
from pybedtools import BedTool

import tools.regions as regtools
import utils
from analysis.motifs import Motifs
from preflights import check_reference_format

STANDARD_GENOMES = ['hg19', 'b37', 'GRCh38', 'mm10']

REFERENCE_VERSION = {
    'hg19': '1.2',
    'b37': '1.2',
    'GRCh38': '1.2',
    'mm10': '1.2',
}

FASTA_SOURCE = {
    'hg19': 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
    'b37': 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
    'GRCh38': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz',
    'mm10': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.5_GRCm38.p3/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz',
}

ASSEMBLY_ACCESSION = {
    'hg19': 'GCA_000001405.1',
    'b37': 'GCA_000001405.1',
    'GRCh38': 'GCA_000001405.15',
    'mm10': 'GCA_000001635.5',
}

GTF_SOURCE = {
    'hg19': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.basic.annotation.gtf.gz',
    'b37': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.basic.annotation.gtf.gz',
    'GRCh38': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.basic.annotation.gtf.gz',
    'mm10': 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.basic.annotation.gtf.gz'
}

ANNOTATION_VERSION = {
    'hg19': 'gencode.v28lift37.basic',
    'b37': 'gencode.v28lift37.basic',
    'GRCh38': 'gencode.v28.basic',
    'mm10': 'gencode.vM17.basic'
}

DHS_SOURCE = {
    'hg19': 'https://www.encodeproject.org/files/ENCFF257KON/@@download/ENCFF257KON.bed.gz',
    'b37': 'https://www.encodeproject.org/files/ENCFF257KON/@@download/ENCFF257KON.bed.gz',
    'GRCh38': 'http://cf.10xgenomics.com/supp/cell-atac/dnase.bed',  # host kundaje's dnase file in 10x website
    'mm10': 'https://www.encodeproject.org/files/ENCFF203SGZ/@@download/ENCFF203SGZ.bed.gz'
}

BLACKLIST_SOURCE = {
    'hg19': 'https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz',
    'b37': 'https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz',
    'GRCh38': 'https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz',
    'mm10': 'https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz'
}


REG_ENSEMBL = {
    'hg19': 'ftp://ftp.ensembl.org/pub/grch37/release-95/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz',
    'b37': 'ftp://ftp.ensembl.org/pub/grch37/release-95/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz',
    'GRCh38': 'ftp://ftp.ensembl.org/pub/release-95/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190122.gff.gz',
    'mm10': 'ftp://ftp.ensembl.org/pub/release-95/regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz'
}

ORGANISM = {
    'hg19': 'Homo_sapiens',
    'b37': 'Homo_sapiens',
    'GRCh38': 'Homo_sapiens',
    'mm10': 'Mus_musculus'
}

# hg19 genome fasta does not have "chr" prefix in major contigs, based on the source file
CONTIG_DEFS = {
    'hg19': {
        'primary_contigs': [str(i) for i in range(1, 23)] + ['X', 'Y'],
        'sex_chromosomes': ['X', 'Y'],
        'non_nuclear_contigs': ['M'],
        'prefix': 'chr'
    },
    'b37': {
        'primary_contigs': [str(i) for i in range(1, 23)] + ['X', 'Y'],
        'sex_chromosomes': ['X', 'Y'],
        'non_nuclear_contigs': ['MT'],
        'prefix': ''
    },
    'GRCh38': {
        'primary_contigs': [str(i) for i in range(1, 23)] + ['X', 'Y'],
        'sex_chromosomes': ['X', 'Y'],
        'non_nuclear_contigs': ['M'],
        'prefix': 'chr',
    },
    'mm10': {
        'primary_contigs': [str(i) for i in range(1, 20)] + ['X', 'Y'],
        'sex_chromosomes': ['X', 'Y'],
        'non_nuclear_contigs': ['M'],
        'prefix': 'chr'
    },
}

# motif collection is from JASPAR 2018 core non-redundant
MOTIF_COLLECTION = {
    'hg19': 'http://cf.10xgenomics.com/supp/cell-atac/jaspar2018_core_vertebrates_non-redundant.pfm',
    'b37': 'http://cf.10xgenomics.com/supp/cell-atac/jaspar2018_core_vertebrates_non-redundant.pfm',
    'GRCh38': 'http://cf.10xgenomics.com/supp/cell-atac/jaspar2018_core_vertebrates_non-redundant.pfm',
    'mm10': 'http://cf.10xgenomics.com/supp/cell-atac/jaspar2018_core_vertebrates_non-redundant.pfm'
}

STANDARD_REF_FILES = {
    'blacklist.bed',
    'ctcf.bed',
    'dnase.bed',
    'enhancer.bed',
    'genes.gtf',
    'genome.fa',
    'promoter.bed',
    'regulatory.gff',
    'tss.bed'
}

class GxfParser:

    def __init__(self, ref_path, genes_in, contig_defs):

        self.ref_path = ref_path
        self.genes_in = genes_in
        self.contig_defs = contig_defs
        self.genes_out = None

    def get_gtf(self):

        self.genes_out = os.path.join(self.ref_path, "genes", "genes.gtf")
        utils.fetch_from_source(self.genes_in, self.genes_out, 'gene annotation')

        # check if genes_out is in GTF format. If in GFF3, convert it to GTF
        if regtools.gtf_or_gff3(self.genes_out) == 'GFF':
            genes_out_gtf = self.genes_out + '.gtf'
            regtools.convert_gff3_to_gtf(self.genes_out, genes_out_gtf)
            shutil.move(genes_out_gtf, self.genes_out)

        # clean chromosome names and filter out non-primary contigs (chrM is not filtered out)
        if self.contig_defs['reference_name'] in STANDARD_GENOMES:
            regtools.clean_chr_name_file(self.genes_out, self.contig_defs, self.ref_path)
        else:
            # check gtf contig names
            regtools.validate_contig_names(self.genes_out, self.contig_defs)

    def get_tss(self):
        # output tss and transcripts bed
        print("Writing TSS and transcripts bed file...")
        sys.stdout.flush()

        tss_out = os.path.join(self.ref_path, 'regions', 'tss.bed')
        tx_out = os.path.join(self.ref_path, 'regions', 'transcripts.bed')
        faidx_file = os.path.join(self.ref_path, 'fasta', 'genome.fa.fai')

        if os.path.exists(tss_out):
            print("{} file already exists! Skip parsing from gene annotation...\n".format(tss_out))
            # When tss.bed exists, check format. If pass, skip tss the generation step. same for transcripts.bed
            regtools.bed_format_checker(tss_out, faidx_file)

            count = utils.quick_line_count(tss_out)
            print("    Parsed {} unique TSS.".format(count))

            if os.path.exists(tx_out):
                regtools.bed_format_checker(tx_out, faidx_file)

                count = utils.quick_line_count(tx_out)
                print("    Parsed {} unique transcripts.".format(count))
            print("done\n")
            return

        with open(tss_out, 'w') as tss_writer, open(tx_out, 'w') as tx_writer:

            has_tx = False
            valid_contig_names = utils.extact_contig_names_from_def(self.ref_path)

            gtf_iter = regtools.bed_format_chcker_iter(self.genes_out)

            for i, fields in enumerate(gtf_iter):

                if self._is_valid_transcript_record(fields, valid_contig_names):
                    has_tx = True

                    # extract TSS and transcript coordinates
                    if fields['strand'] == '+':
                        tss_start = fields['start']
                        tss_end = fields['start'] + 1

                        tx_start = fields['start']
                        tx_end = fields['end']

                    elif fields['strand'] == '-':
                        tss_start = fields['end'] - 1
                        tss_end = fields['end']

                        tx_start = fields['start']
                        tx_end = fields['end']

                    # check gene_name is one of the keys in fields
                    try:
                        name = fields['gene_name']
                    except AttributeError:
                        msg = "Invalid GTF entry at line {} of gene annotation: cannot find \"gene_name\" in the attribute column.".format(i+1)
                        sys.exit(msg)

                    if 'gene_type' in fields.attrs:
                        tss_row = [fields['chrom'], str(tss_start), str(tss_end), name, '.', fields['strand'], fields['gene_type']]
                        tx_row = [fields['chrom'], str(tx_start), str(tx_end), name, '.', fields['strand'], fields['gene_type']]
                    else:
                        tss_row = [fields['chrom'], str(tss_start), str(tss_end), name, '.', fields['strand']]
                        tx_row = [fields['chrom'], str(tx_start), str(tx_end), name, '.', fields['strand']]

                    tss_writer.write('\t'.join(tss_row) + '\n')
                    tx_writer.write('\t'.join(tx_row) + '\n')

        if not has_tx:
            sys.exit("Invalid gene annotation input: "
                     "cannot find \"transcript\" or  \"mRNA\" type in the 3rd column.")

        # print("done parsing tss..."
        regtools.sort_and_uniq_bed(tss_out, self.ref_path)
        regtools.sort_and_uniq_bed(tx_out, self.ref_path)

        regtools.bed_format_checker(tss_out, faidx_file)
        regtools.bed_format_checker(tx_out, faidx_file)

        count_tss = utils.quick_line_count(tss_out)
        count_tx = utils.quick_line_count(tx_out)

        print("    Parsed {} unique TSS and {} unique transcripts.\ndone\n".format(count_tss, count_tx))

    @staticmethod
    def _is_valid_transcript_record(bed_fields, valid_contig_names):
        """test the row of the bed file is a valid record for a transcript"""

        cond1_is_tx = bed_fields[2] == "transcript" or bed_fields[2] == "mRNA"
        cond2_valid_contig = bed_fields['chrom'] in valid_contig_names

        # filter PAR genes for human genome
        # when transcript_id is not present, skip the filtering.
        # covers standard genome, human gene annotation source from
        # GENCODE and ensembl
        cond3_not_par = not bed_fields.attrs["transcript_id"].endswith(
            "_PAR_Y") if "transcript_id" in bed_fields.attrs else True

        return all([cond1_is_tx, cond2_valid_contig, cond3_not_par])


class EnsemblREG(GxfParser):

    def get_reg_elements(self):

        # fetch regulatory build gff from source
        self.genes_out = self.get_gff()

        # parse gff to generate regulatory regions bed
        promoter_out = os.path.join(self.ref_path, 'regions', 'promoter.bed')
        enhancer_out = os.path.join(self.ref_path, 'regions', 'enhancer.bed')
        ctcf_out = os.path.join(self.ref_path, 'regions', 'ctcf_raw.bed')

        print("Parsing functional annotation bed files...")
        count = {'promoter': 0,
                 'enhancer': 0,
                 'ctcf': 0}

        with open(promoter_out, 'w') as pro, open(enhancer_out, 'w') as enh, open(ctcf_out, 'w') as ctcf:

            promoter_writer = csv.writer(pro, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            enhancer_writer = csv.writer(enh, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            ctcf_writer = csv.writer(ctcf, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')

            gtf_iter = regtools.bed_format_chcker_iter(self.genes_out)
            for fields in gtf_iter:

                # 1based coord. from gtf/gff is accommodated in reader_iter
                reg_start = fields['start']
                reg_end = fields['end']
                row = [fields['chrom'], reg_start, reg_end, '.', '.', fields['strand']]

                if fields['feature_type'] == 'Promoter':
                    count['promoter'] += 1
                    promoter_writer.writerow(row)

                elif fields['feature_type'] == 'CTCF Binding Site':
                    count['ctcf'] += 1
                    ctcf_writer.writerow(row)

                else:
                    count['enhancer'] += 1
                    enhancer_writer.writerow(row)

        regtools.sort_and_uniq_bed(promoter_out, self.ref_path)
        regtools.sort_and_uniq_bed(enhancer_out, self.ref_path)
        regtools.sort_and_uniq_bed(ctcf_out, self.ref_path)

        for value, key in count.items():
            print("    Parsed {} {} regions.".format(key, value))
        print("done")

        # scan CTCF motifs in CTCF sites
        print("\nRefine CTCF sites...")
        ctcf_count = self.refine_ctcf(ctcf_out)
        print("    Refined CTCF sites to {} sites...\ndone\n".format(ctcf_count))

    def get_gff(self):

        regulatory_out = os.path.join(self.ref_path, "genes", "regulatory.gff")
        utils.fetch_from_source(self.genes_in, regulatory_out, 'regulatory element annotation')

        # clean chromosome names and filter out non-primary contigs (chrM is not filtered out)
        regtools.clean_chr_name_file(regulatory_out, self.contig_defs, self.ref_path)

        # sorted gff with numerical chromosome order.
        regtools.sort_and_uniq_bed(regulatory_out, self.ref_path)

        # validate formatting
        faidx_file = os.path.join(self.ref_path, 'fasta', 'genome.fa.fai')
        regtools.bed_format_checker(regulatory_out, faidx_file)

        return regulatory_out

    def refine_ctcf(self, ctcf_out):
        motifs = Motifs(self.ref_path)
        ctcf_refined = os.path.join(self.ref_path, 'regions', 'ctcf.bed')

        print("Scan CTCF motif in CTCF sites...")
        ctcf_bedtool = BedTool(ctcf_out)
        motifs.scan_motif_from_bed(ctcf_bedtool, 'CTCF_MA0139.1', ctcf_refined, out_format="motif", use_genome_bg=True, pseudocount=1, pvalue=1e-5, window=7)
        regtools.sort_and_uniq_bed(ctcf_refined, self.ref_path)
        subprocess.check_call(['rm', ctcf_out])  # delete the raw sites file

        count = utils.quick_line_count(ctcf_refined)
        return count


class FastaBuilder:
    def __init__(self, ref_path, fasta_in, contig_defs):
        self.ref_path = ref_path
        self.fasta_in = fasta_in
        self.fasta_out = None
        self.contig_defs = contig_defs

    # TODO: when user supplies custom fasta file, need to create a fasta validator

    def get_fasta(self):

        self.fasta_out = os.path.join(self.ref_path, "fasta", "genome.fa")
        utils.fetch_from_source(self.fasta_in, self.fasta_out, 'fasta')

    def rename_contig(self):
        """
        rename genome contig names: adding "chr" prefix, only used when building hg19
        :return: a new genome.fa file with new contig names
        """

        fasta_renamed = self.fasta_out + '.tmp'

        with open(self.fasta_out) as f, open(fasta_renamed, 'w') as renamed:
            new_chr_name = None
            for line in f:
                if line.startswith('>'):
                    chr_name = line.strip().lstrip('>').split()[0]
                    new_chr_name = regtools.clean_chr_name(chr_name, self.contig_defs, keep_nonstd=True)
                    if new_chr_name is not None:
                        renamed.write('>' + new_chr_name + '\n')
                else:
                    if new_chr_name is not None:
                        renamed.write(line)

        shutil.move(fasta_renamed, self.fasta_out)

    def get_contig_defs(self, species_prefixes=[""]):
        contig_defs = self.contig_defs

        if contig_defs['reference_name'] in STANDARD_GENOMES:

            if contig_defs['sex_chromosomes'] == ['X', 'Y']:

                sex_chromosomes = {
                    "_male": {
                        contig_defs['prefix'] + "X": 1,
                        contig_defs['prefix'] + "Y": 1
                    },
                    "_female": {
                        contig_defs['prefix'] + "X": 2,
                        contig_defs['prefix'] + "Y": 0
                    }
                }
            else:
                sex_chromosomes = contig_defs['sex_chromosomes']

            contig_defs_to_json = [
                ('species_prefixes', species_prefixes),
                ('primary_contigs', [contig_defs['prefix'] + i for i in contig_defs['primary_contigs']]),
                ('sex_chromosomes', sex_chromosomes),
                ("non_nuclear_contigs", [contig_defs['prefix'] + i for i in contig_defs['non_nuclear_contigs']])
            ]
        else:
            contig_defs_to_json = [
                ('species_prefixes', species_prefixes),
                ('primary_contigs', contig_defs['PRIMARY_CONTIGS']),
                ('non_nuclear_contigs', contig_defs['NON_NUCLEAR_CONTIGS'])
            ]
            for key in ['sex_chromosomes']:
                if key in contig_defs:
                    contig_defs_to_json.append((key, contig_defs[key]))

        contig_defs_to_json = collections.OrderedDict(contig_defs_to_json)
        contig_defs_out = os.path.join(self.ref_path, 'fasta', 'contig-defs.json')
        with open(contig_defs_out, 'w') as f:
            f.write(json.dumps(contig_defs_to_json, indent=4))

    def generate_index(self):
        print("Generating samtools index...")
        sys.stdout.flush()
        try:
            subprocess.check_call(['samtools', 'faidx', self.fasta_out])
        except subprocess.CalledProcessError:
            sys.exit("Invalid genome fasta input, please check the source.")

        print("done\n")

        # for custom references, validate fasta contig names match definition in contig_defs:
        # PRIMARY_CONTIGS must be a subset of contig names in genome.fa
        if self.contig_defs['reference_name'] not in STANDARD_GENOMES:
            regtools.validate_contig_names(self.fasta_out + '.fai', self.contig_defs)

        print("Generating pyfasta indexes...")
        sys.stdout.flush()
        pyf = pyfasta.Fasta(self.fasta_out, key_fn=lambda x: x.split()[0])
        contigs = len(pyf)
        size = sum(len(pyf[contig]) for contig in pyf)

        print("    Number of contigs: %d\n    Total genome size: %d" % (contigs, size))
        print("done\n")

    def generate_bwa_index(self):

        print("Generating bwa index (may take over an hour for a 3Gb genome)...")
        sys.stdout.flush()
        subprocess.check_call(['bwa', 'index', self.fasta_out])
        print("done\n")


class ReferenceBuilder():
    def __init__(self, ref_name, config=None):
        """
        Initiation of ReferenceBuilder
        :param ref_name: name of the genome
        :param config: the path of the configuration file. Ignored with building standard genomes
        """

        self.ref_name = ref_name
        self.ref_path = None

        if self.ref_name in STANDARD_GENOMES:
            self.config = self.parse_ref_config(CONTIG_DEFS[self.ref_name])
            self.fasta_source = FASTA_SOURCE[self.ref_name]
            self.genes_source = GTF_SOURCE[self.ref_name]
            self.regulatory_source = REG_ENSEMBL[self.ref_name]

        else:
            self.config = self.parse_ref_config(config)
            self.fasta_source = self.config['GENOME_FASTA_INPUT']
            self.genes_source = self.config['GENE_ANNOTATION_INPUT']
            self.regulatory_source = None

    def build_reference(self):

        # Step 0: folder structure
        self.ref_path = self.make_ref_dir(self.ref_name)

        # Step 1: Write a metadata.json file to the reference directory with Species, assembly, and annotation metadata.
        self.create_metadata()

        # Step 2: Genome sequences and indexes
        genome_seq = FastaBuilder(self.ref_path, self.fasta_source, self.config)

        # Download fasta files
        genome_seq.get_fasta()

        if self.ref_name == 'hg19':
            genome_seq.rename_contig()

        # generate contig-defs.json
        genome_seq.get_contig_defs()

        # generate indexes of genome sequences
        genome_seq.generate_index()

        # Step 3 gene annotation GTF and TSS
        # Check the existence of gtf files, if not, download gtf files
        anno_gtf = GxfParser(self.ref_path, self.genes_source, self.config)
        anno_gtf.get_gtf()

        # generate TSS bed file
        anno_gtf.get_tss()

        # do the bwa index last
        genome_seq.generate_bwa_index()

        # Step 4: Regulatory element (Promoter, Enhancer, CTCF sites, DHS and Blacklist)
        # Download motif database
        self.get_motif()

        if self.ref_name in STANDARD_GENOMES:
            # Promoter, Enhancer and CTCF regions
            anno_reg = EnsemblREG(self.ref_path, self.regulatory_source, self.config)
            anno_reg.get_reg_elements()

            # DNase hypersensitivity sites
            self.get_DHS()

            # Black list
            self.get_blacklist()

        self.clean_permissions(self.ref_path)

        # Step 5: check the newly built reference
        check_reference_format(self.ref_path)

    @staticmethod
    def make_ref_dir(ref_name):
        # Create the target reference folder path
        ref_folder = ref_name
        ref_path = os.path.join(os.getcwd(), ref_folder)

        # Check for write permissions in current directory
        if not os.access(os.getcwd(), os.W_OK):
            sys.exit("You do not have write permission in the current working directory.")

        # Check that destination folder doesn't already exist
        if os.path.exists(ref_path):
            print("Destination reference folder already exists: {}\n".format(ref_path))

        else:
            # Create reference folder
            print("Creating new reference folder at {}".format(ref_path))
            os.mkdir(ref_path)

            # Write out genome identifier
            with open(os.path.join(ref_path, "genome"), "w") as f:
                f.write(ref_name + "\n")

            # Create reference folder structure
            os.mkdir(os.path.join(ref_path, "fasta"))
            os.mkdir(os.path.join(ref_path, "genes"))
            os.mkdir(os.path.join(ref_path, "regions"))

        return ref_path

    def parse_ref_config(self, config_in):
        """

        :param config_in: a pre-defined dictionary if building a standard reference,
                          or a file path for configuration if building a custom reference
        :return: a dictionary of contig definition (for standard reference) plus metadata (for custom reference)
        """

        if self.ref_name in STANDARD_GENOMES:
            ref_config = config_in
            ref_config['reference_name'] = self.ref_name

        else:
            if config_in is None:
                sys.exit("Error: reference configuration file is required for custom references.")

            with open(config_in) as cfg:
                config_text = cfg.read()

            try:
                ref_config = hjson.loads(config_text)
            except hjson.scanner.HjsonDecodeError as e:
                note = "\nNote: do not use the escape character '\\' in input file paths."
                sys.exit(
                    "Error in parsing of the reference configuration file: " + str(e) + note
                )

            ref_config['reference_name'] = self.ref_name

            # validate input
            try:
                ref_config['GENOME_FASTA_INPUT'] = str(ref_config['GENOME_FASTA_INPUT'])
                ref_config['GENE_ANNOTATION_INPUT'] = str(ref_config['GENE_ANNOTATION_INPUT'])
                ref_config['MOTIF_INPUT'] = str(ref_config['MOTIF_INPUT'])
                ref_config['ORGANISM'] = str(ref_config['ORGANISM'])
                ref_config['PRIMARY_CONTIGS'] = [str(i) for i in ref_config['PRIMARY_CONTIGS']]

                assert ref_config['PRIMARY_CONTIGS'] != [] and ref_config['PRIMARY_CONTIGS'] != [""], "PRIMARY_CONTIGS can not be empty."

                ref_config['NON_NUCLEAR_CONTIGS'] = [str(i) for i in ref_config['NON_NUCLEAR_CONTIGS']]
                if ref_config['NON_NUCLEAR_CONTIGS'] == [""] or ref_config['NON_NUCLEAR_CONTIGS'] is None:
                    ref_config['NON_NUCLEAR_CONTIGS'] = []

            except:
                sys.exit("Error: custom reference configuration file is invalid.")

        return ref_config

    def create_metadata(self):
        metadata = {}
        if self.ref_name in STANDARD_GENOMES:
            metadata["assembly"] = self.ref_name
            metadata["organism"] = ORGANISM[self.ref_name]
            metadata["assembly_accession"] = ASSEMBLY_ACCESSION[self.ref_name]
            metadata["assembly_fasta_url"] = FASTA_SOURCE[self.ref_name]
            metadata["annotation"] = ANNOTATION_VERSION[self.ref_name]
            metadata["annotation_gtf_url"] = GTF_SOURCE[self.ref_name]
            metadata["version"] = REFERENCE_VERSION[self.ref_name]
        else:
            metadata["assembly"] = "custom"
            metadata["organism"] = self.config['ORGANISM']
            metadata["assembly_accession"] = "custom"
            metadata["assembly_fasta_url"] = self.config['GENOME_FASTA_INPUT']
            metadata["annotation"] = "custom"
            metadata["annotation_gtf_url"] = self.config['GENE_ANNOTATION_INPUT']
            metadata["version"] = "custom"

        with open(os.path.join(self.ref_path, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=4)

    def get_DHS(self):
        """
        download and clean the dnase.bed. Only called when creating standard genomes.
        :return:
        """

        source = DHS_SOURCE[self.ref_name]
        path = os.path.join(self.ref_path, 'regions', 'dnase.bed')
        utils.download_from_url(source, path)

        regtools.clean_chr_name_file(path, self.config, self.ref_path)
        regtools.sort_and_uniq_bed(path, self.ref_path)

    def get_blacklist(self):

        source = BLACKLIST_SOURCE[self.ref_name]
        path = os.path.join(self.ref_path, 'regions', 'blacklist.bed')
        utils.download_from_url(source, path)

        # hg19 blacklist.bed has a coordinate overflow issue for the chrM row.
        # The entire chrM is set to be blacklist but the annotation does not match the contig length
        regtools.clean_chr_name_file(path, self.config, self.ref_path)
        regtools.sort_and_uniq_bed(path, self.ref_path)

    def get_motif(self):
        if self.ref_name in STANDARD_GENOMES:
            source = MOTIF_COLLECTION[self.ref_name]
        elif self.config['MOTIF_INPUT'] is not None and self.config['MOTIF_INPUT'] != "":
            source = self.config['MOTIF_INPUT']
        else:
            return

        path = os.path.join(self.ref_path, 'regions', 'motifs.pfm')
        utils.fetch_from_source(source, path, 'pfm')

        # clean motif names to remove tabs
        self.clean_motif_names(path)

    # clean motif names
    def clean_motif_names(self, motif_infile):
        """
        Clean motif naming scheme and check motif file formatting
        :param motif_infile: file path to the motif file
        :return: a valid motif file. if there's issue in the formatting, sys.exit() with an error message
        """

        try:
            # for custom ref, motif.pfm can be formatted as the final form
            utils.motif_format_checker(motif_infile)
            return

        except SystemExit:
            motif_in = open(motif_infile, 'r')
            motif_outfile = os.path.join(self.ref_path, 'regions', motif_infile + '.tmp')

            with open(motif_outfile, 'w') as out:
                for i, line in enumerate(motif_in):
                    if line.startswith('>'):
                        motif_name_fields = line.strip().lstrip('>').split()
                        if motif_name_fields == [""]:
                            sys.exit("Invalid motif file: missing motif names at line {}".format(i+1))
                        else:
                            if len(motif_name_fields) == 1:
                                print("Warning: No gene name is detected in the motif name at line {}.".format(i+1))
                                print("         If {} contains a motif id and a gene name, they should be separated by a tab (\\t).".format(motif_name_fields[0]))
                                motif_name = ">{}\n".format(motif_name_fields[0])
                            else:
                                motif_version = motif_name_fields[0]
                                motif_gene = motif_name_fields[1]
                                motif_name = ">{}_{}\n".format(motif_gene, motif_version)
                        out.write(motif_name)
                    else:
                        out.write(line)
            shutil.move(motif_outfile, motif_infile)

            utils.motif_format_checker(motif_infile)

    # clean permissions
    @staticmethod
    def clean_permissions(ref_path):
        import stat
        print("Finishing up...")

        for root, subdirs, files in os.walk(ref_path):

            perm = stat.S_IRUSR + stat.S_IWUSR + stat.S_IRGRP + stat.S_IWGRP + stat.S_IROTH
            root = os.path.abspath(root)
            for fname in files:
                f_path = os.path.join(root, fname)
                os.chmod(f_path, perm)

                # touch genome.fa.* for the pyfasta indexing issue
                if fname.startswith('genome.fa.'):
                    os.utime(f_path, None)


class BarnyardBuilder():
    def __init__(self, human_build, mouse_build):
        self.human_build = human_build
        self.mouse_build = mouse_build

        self.barnyard_name = human_build + '_and_' + mouse_build
        self.barnyard_path = None

    def build_reference(self):

        # Step 1: folder structure
        self.barnyard_path = ReferenceBuilder.make_ref_dir(self.barnyard_name)

        # Step 1: Write a metadata.json file to the reference directory with Species, assembly, and annotation metadata.
        metadata = self.combine_metadata()
        with open(os.path.join(self.barnyard_path, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=4)

        # Step 2: extract referece files from human and mouse reference, which are built previously
        human_ref_path = os.path.join(os.getcwd(), self.human_build)
        mouse_ref_path = os.path.join(os.getcwd(), self.mouse_build)

        human_ref_files = self.get_ref_files(human_ref_path)
        mouse_ref_files = self.get_ref_files(mouse_ref_path)

        # Step 3, concatenate match files and update chromosome names
        for h, m in zip(human_ref_files, mouse_ref_files):
            assert os.path.basename(h) == os.path.basename(m)

            target_file = h.split('/')
            target_file = os.path.join(self.barnyard_path, target_file[-2], target_file[-1])

            if h.endswith('.fa'):
                self.merge_from_two_sources(h, m, target_file, 'fasta')
            else:
                self.merge_from_two_sources(h, m, target_file, 'tab')

        # Step 4, build contig-defs.json
        h_contigs = CONTIG_DEFS[self.human_build]
        m_contigs = CONTIG_DEFS[self.mouse_build]
        h_prefix = self.human_build + '_' + h_contigs['prefix']
        m_prefix = self.mouse_build + '_' + m_contigs['prefix']

        contig_defs = {
            'reference_name': self.barnyard_name,
            'PRIMARY_CONTIGS': [h_prefix + i for i in h_contigs['primary_contigs']] +
                               [m_prefix + j for j in m_contigs['primary_contigs']],
            'sex_chromosomes': {
                self.human_build + "_male": {
                    h_prefix + "X": 1,
                    h_prefix + "Y": 1
                },
                self.human_build + "_female": {
                    h_prefix + "X": 2,
                    h_prefix + "Y": 0
                },
                self.mouse_build + "_male": {
                    m_prefix + "X": 1,
                    m_prefix + "Y": 1
                },
                self.mouse_build + "_female": {
                    m_prefix + "X": 2,
                    m_prefix + "Y": 0
                }
            },
            'NON_NUCLEAR_CONTIGS': [h_prefix + i for i in h_contigs['non_nuclear_contigs']] +
                                   [m_prefix + i for i in m_contigs['non_nuclear_contigs']],
            'prefix': ''
        }

        fa_file = os.path.join(self.barnyard_path, 'fasta', 'genome.fa')
        genome_seq = FastaBuilder(self.barnyard_path, fa_file, contig_defs=contig_defs)
        genome_seq.get_fasta()
        genome_seq.get_contig_defs(species_prefixes=[self.human_build, self.mouse_build])

        # Step 5, build sequence indexes
        genome_seq.generate_index()
        genome_seq.generate_bwa_index()

        # Step 6, clean permissions
        ReferenceBuilder.clean_permissions(self.barnyard_path)

    def combine_metadata(self):

        human_metadata = json.load(open(os.path.join(self.human_build, "metadata.json")))
        mouse_metadata = json.load(open(os.path.join(self.mouse_build, "metadata.json")))

        # Step 1: Write a metadata.json file to the reference directory with Species, assembly, and annotation metadata.
        metadata = {}
        metadata["assembly"] = human_metadata["assembly"] + "_and_" + mouse_metadata["assembly"]
        metadata["organism"] = human_metadata["organism"] + "_and_" + mouse_metadata["organism"]

        metadata["assembly_accession"] = "barnyard"
        metadata["assembly_fasta_url"] = "barnyard"
        metadata["annotation"] = "barnyard"
        metadata["annotation_gtf_url"] = "barnyard"
        metadata["version"] = human_metadata["version"] + "_and_" + mouse_metadata["version"]
        return metadata

    def get_ref_files(self, ref_path):

        print("Collect reference data from {}...\n".format(ref_path))

        path_iter = os.walk(ref_path)
        path_iter.next()
        saved = []
        seen = set()
        for parent, folders, files in path_iter:
            for file in files:
                if file in STANDARD_REF_FILES:
                    saved.append(os.path.join(parent, file))
                    seen.add(file)

        if not seen == STANDARD_REF_FILES:
            sys.exit('Reference data at {} is not built properly, missing files of {}'.format(ref_path, list(STANDARD_REF_FILES - seen)))

        saved.sort()
        return saved

    def merge_from_two_sources(self, source1, source2, file_out, file_type):
        # source1: human source
        # source2: mouse source
        # file_out: file name of the output file
        # file_type: fasta or tab

        prefix1 = self.human_build
        prefix2 = self.mouse_build

        print('Generating {} from {} and {}...'.format(file_out, prefix1, prefix2))

        with open(file_out, 'w') as f:
            if file_type == 'fasta':
                self._fasta_rename(source1, f, prefix1)
                self._fasta_rename(source2, f, prefix2)

            if file_type == 'tab':
                self._tab_rename(source1, f, prefix1)
                self._tab_rename(source2, f, prefix2)

        print("done\n")

    @staticmethod
    def _fasta_rename(fasta_in, fasta_out, prefix):
        # input a fasta file, add prefix to chromosome identifiers
        with open(fasta_in) as f:
            for line in f:
                if line.startswith(">"):
                    fields = line[1:].strip().split()
                    new_name = '>' + prefix + '_' + fields[0] + '\n'
                    fasta_out.write(new_name)
                else:
                    fasta_out.write(line)

    @staticmethod
    def _tab_rename(file_in, file_out, prefix):
        # input a tab-delimited file in which the first column is the chromosome, add prefix to it
        with open(file_in) as f:
            for line in f:
                fields = line.split('\t')
                fields[0] = prefix + '_' + fields[0]
                file_out.write("\t".join(fields))
