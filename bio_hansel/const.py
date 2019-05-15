# -*- coding: utf-8 -*-

import re
from pkg_resources import resource_filename

from bio_hansel import program_name
from bio_hansel.subtyping_params import SubtypingParams

SCHEME_FASTAS = {'heidelberg': {'file': resource_filename(program_name, 'data/heidelberg/kmers.fasta'),
                                'version': '0.5.0',
                                'subtyping_params': SubtypingParams(low_coverage_depth_freq=20)},
                 'enteritidis': {'file': resource_filename(program_name, 'data/enteritidis/kmers.fasta'),
                                 'version': '1.0.7',
                                 'subtyping_params': SubtypingParams(low_coverage_depth_freq=50)},
                 'typhi': {'file': resource_filename(program_name, 'data/typhi/kmers.fasta'),
                                 'version': '1.1.0',
                                 'subtyping_params': SubtypingParams(low_coverage_depth_freq=20)},
                 'tb_speciation': {'file': resource_filename(program_name, 'data/m.tuberculosis/kmers.fasta'),
                                 'version': '1.0.1;',
                                 'subtyping_params': SubtypingParams(low_coverage_depth_freq=20)},
                 'typhimurium': {'file': resource_filename(program_name, 'data/typhimurium/kmers.fasta'),
                                 'version': '0.5.5;',
                                 'subtyping_params': SubtypingParams(low_coverage_depth_freq=20)}}


bases_dict = {
'A': ['A'],
'C': ['C'],
'G': ['G'],
'T': ['T'],
'R': ['A', 'G'],
'Y': ['C', 'T'],
'S': ['G', 'C'],
'W': ['A', 'T'],
'K': ['G', 'T'],
'M': ['A', 'C'],
'B': ['C', 'G', 'T'],
'D': ['A', 'G', 'T'],
'H': ['A', 'C', 'T'],
'V': ['A', 'C', 'G'],
'N': ['A', 'C', 'G', 'T'],}


COLUMNS_TO_REMOVE = '''
pident
length
mismatch
gapopen
qstart
qend
sstart
send
evalue
bitscore
qlen
slen
sseq
is_trunc
coverage
'''.strip().split('\n')

# These are present within the subtype module.
SUBTYPE_SUMMARY_COLS = """
sample
scheme
scheme_version
subtype
all_subtypes
kmers_matching_subtype
are_subtypes_consistent
inconsistent_subtypes
n_kmers_matching_all
n_kmers_matching_all_expected
n_kmers_matching_positive
n_kmers_matching_positive_expected
n_kmers_matching_subtype
n_kmers_matching_subtype_expected
file_path
avg_kmer_coverage
qc_status
qc_message
""".strip().split('\n')


SIMPLE_SUMMARY_COLS = """
sample
subtype
coverage
qc_status
qc_message
""".strip().split('\n')

REGEX_FASTQ = re.compile(r'^(.+)\.(fastq|fq|fastqsanger)(\.gz)?$')
REGEX_FASTA = re.compile(r'^.+\.(fasta|fa|fna|fas)(\.gz)?$')

JSON_EXT_TMPL = '{}.json'
