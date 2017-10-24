# -*- coding: utf-8 -*-

FASTA_COLUMNS_TO_REMOVE = '''
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
seq
coverage
is_trunc
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
confident_is_subtype
reached_min_tiles""".strip().split('\n')

SIMPLE_SUMMARY_COLS = """
sample
subtype
result
""".strip().split('\n')

MIXED_SUBTYPE_ERROR = "ERROR: Mixed Subtypes"
OK_SUBTYPE = "PASS: Subtypes not mixed"
INSUFFICIENT_NUM_TILES = "ERROR: Insufficient number of SNV targets found!"
OK_NUM_TILES = "PASS: Expected number of tiles reached"
