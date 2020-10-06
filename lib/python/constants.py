"""
These are constants that are specific to this pipeline.  Other
general constants can be found in tenkit.constants.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division

# Product name, used in bam header lines
TENX_PRODUCT_NAME = 'cellranger-atac'

# Sequence names used in bcl2fastq/mkfastq output files
BCL2FASTQ_SEQNAMES = {'read1': 'R1', 'read2': 'R3',
                      'barcode': 'R2', 'sample_index': 'I1',
                      'R1': 'R1', 'R2': 'R3'}

# Threshold for calling a read low mapping quality
LOW_MAPQ_THRESHOLD = 30

# These override the default values from tenkit
RAW_BARCODE_TAG = 'CR'
PROCESSED_BARCODE_TAG = 'CB'
RAW_BARCODE_QUAL_TAG = 'CY'

# Additional tags for the 5' global position
# of read and its mate. Used for sorting and
# marking duplicates
SELF_FIVE_PRIME_POS_TAG = 'GP'
MATE_FIVE_PRIME_POS_TAG = 'MP'
MATE_MAPPING_QUALITY_TAG = 'MQ'

# Downsampled number of reads to report standardized metrics at
DOWNSAMPLED_READS = 30000000
RPC_50K = 50000
RPC_30K = 30000
RPC_10K = 10000

# Limit the total number of cells called for each species to no more than this
MAXIMUM_CELLS_PER_SPECIES = 20000
MINIMUM_COUNT = 5
WHITELIST_CONTAM_RATE = 0.02

NUM_READS_ALIGN_CONTAMINATION = 10000
NUM_CHUNKS_CONTAMINATION = 100
GENOMES_CONTAMINATION = ['hg19', 'sacCer3', 'ecoli_K12_MG1655', 'ecoli_K12_DH10B', 'PhiX', 'staph_PCI_1200', 'bac_vectors']

PWM_MATCH_PVAL_THRESHOLD = 1E-7

# some analysis constants
ALLOWED_FACTORIZATIONS = set(['plsa', 'pca', 'lsa'])
CLUSTER_FILE_HEAD = {'pca': 'pca_kmeans', 'lsa': 'lsa_kmeans', 'plsa': 'plsa_kmeans'}
DEFAULT_FACTORIZATION = "lsa"

# aggr
ALLOWED_NORMALIZATIONS = set([None, 'unique', 'signal_mean'])

# CIGAR codes
CIGAR_SOFTCLIP = 4

# Used for position tagging
MAX_I32 = 2147483647

# Ordered valid DNA bases.  Don't change the order!
VALID_BASES = ['A', 'C', 'G', 'T']

# Peak calling parameters
WINDOW_SIZE = 400
# Peaks closer than this will be merged into a single peak.
PEAK_MERGE_DISTANCE = 500
PEAK_ODDS_RATIO = 1 / 5
MAX_INSERT_SIZE = 1200

NO_BARCODE = 'NO_BARCODE'

# for preflight checks on fragments
FRAGMENTS_SCAN_SIZE = 100000

# gene types we care about for annotation
TRANSCRIPT_ANNOTATION_GENE_TYPES = ['protein_coding',
                       'TR_V_gene', 'TR_D_gene', 'TR_J_gene', 'TR_C_gene',
                       'IG_V_gene', 'IG_D_gene', 'IG_J_gene', 'IG_C_gene']

