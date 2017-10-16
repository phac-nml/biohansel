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


MIXED_SUBTYPE_ERROR = "ERROR: MIXED SUBTYPES"
OK_SUBTYPE = "Subtypes not mixed"
INSUFFICIENT_NUM_TILES = "ERROR: Insufficient number of SNV targets found!"
OK_NUM_TILES = "Expected number of tiles reached"
