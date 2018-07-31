# -*- coding: utf-8 -*-

import re
from pkg_resources import resource_filename

from bio_hansel import program_name
from bio_hansel.subtyping_params import SubtypingParams

SCHEME_FASTAS = {'heidelberg': {'file': resource_filename(program_name, 'data/heidelberg/tiles.fasta'),
                                'version': '0.5.0',
                                'subtyping_params': SubtypingParams(low_coverage_depth_freq=20)},
                 'enteritidis': {'file': resource_filename(program_name, 'data/enteritidis/tiles.fasta'),
                                 'version': '0.8.0',
                                 'subtyping_params': SubtypingParams(low_coverage_depth_freq=50)}}

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
tiles_matching_subtype
are_subtypes_consistent
inconsistent_subtypes
n_tiles_matching_all
n_tiles_matching_all_expected
n_tiles_matching_positive
n_tiles_matching_positive_expected
n_tiles_matching_subtype
n_tiles_matching_subtype_expected
file_path
avg_tile_coverage
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
