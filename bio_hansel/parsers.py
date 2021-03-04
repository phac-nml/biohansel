# -*- coding: utf-8 -*-

import logging

import re

VALID_NUCLEOTIDES = {'A', 'a',
                     'C', 'c',
                     'G', 'g',
                     'T', 't',
                     'R', 'r',
                     'Y', 'y',
                     'S', 's',
                     'W', 'w',
                     'K', 'k',
                     'M', 'm',
                     'B', 'b',
                     'D', 'd',
                     'H', 'h',
                     'V', 'v',
                     'N', 'n',
                     'X', 'x', }  # X for masked nucleotides

REGEX_GZIPPED = re.compile(r'^.+\.gz$')


# SimpleFastaParser function from BioPython
# https://github.com/biopython/biopython/blob/92c07ce7dce91078edc9b77fb71dbbe5646565bb/Bio/SeqIO/FastaIO.py#L24
def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.
    Arguments:
     - handle - input stream opened in text mode
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "").upper()
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "").upper()


def parse_fasta(filepath):
    """Parse a FASTA/FASTA.GZ file returning a generator yielding tuples of fasta headers to sequences.

    Args:
        filepath (str): Fasta file path

    Returns:
        generator: yields tuples of (<fasta header>, <fasta sequence>)
    """
    if REGEX_GZIPPED.match(filepath):
        logging.debug('Opening "%s" as gzipped file', filepath)
        # using os.popen with zcat since it is much faster than gzip.open or gzip.open(io.BufferedReader)
        # http://aripollak.com/pythongzipbenchmarks/
        # assumes Linux os with zcat installed
        import os
        with os.popen('zcat < {}'.format(filepath)) as f:
            yield from SimpleFastaParser(f)
    else:
        with open(filepath, 'r') as f:
            yield from SimpleFastaParser(f)


def parse_fastq(filepath):
    """Parse a FASTQ/FASTQ.GZ file returning a generator yielding tuples of FASTQ entry headers and sequences.

    Args:
        filepath (str): FASTQ/FASTQ.GZ file path

    Returns:
        generator: yields tuples of (<fastq header>, <fastq sequence>)
    """
    if REGEX_GZIPPED.match(filepath):
        logging.debug('Opening "%s" as gzipped file', filepath)
        # using os.popen with zcat since it is much faster than gzip.open or gzip.open(io.BufferedReader)
        # http://aripollak.com/pythongzipbenchmarks/
        # assumes Linux os with zcat installed
        import os
        with os.popen('zcat < {}'.format(filepath)) as f:
            yield from _parse_fastq(f)
    else:
        with open(filepath, 'r') as f:
            yield from _parse_fastq(f)


def _parse_fastq(f):
    """Simple FASTQ parser which yields the header and sequence ignoring the quality scores

    Args:
        f: file-like object

    Yields:
        Tuple of FASTQ entry header and sequence
    """
    header = ''
    seq = ''
    skip = False
    for line in f:
        if skip:
            skip = False
            continue
        line = line.strip()
        if line == '':
            continue
        if line[0] == '@':
            header = line.replace('@', '')
        elif line[0] == '+':
            yield header, seq
            skip = True
        else:
            seq = line.upper()
