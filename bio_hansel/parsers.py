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


def parse_fasta(filepath):
    '''Parse a FASTA/FASTA.GZ file returning a generator yielding tuples of fasta headers to sequences.

    Args:
        filepath (str): Fasta file path

    Returns:
        generator: yields tuples of (<fasta header>, <fasta sequence>)
    '''
    if REGEX_GZIPPED.match(filepath):
        logging.debug('Opening "%s" as gzipped file', filepath)
        # using os.popen with zcat since it is much faster than gzip.open or gzip.open(io.BufferedReader)
        # http://aripollak.com/pythongzipbenchmarks/
        # assumes Linux os with zcat installed
        import os
        with os.popen('zcat < {}'.format(filepath)) as f:
            yield from _parse_fasta(f, filepath)
    else:
        with open(filepath, 'r') as f:
            yield from _parse_fasta(f, filepath)


def _parse_fasta(f, filepath):
    seqs = []
    header = ''
    line_count = 0
    for line in f:
        if isinstance(line, bytes):
            line = line.decode()
        line = line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if header == '':
                header = line.replace('>', '')
            else:
                yield header, ''.join(seqs)
                seqs = []
                header = line.replace('>', '')
        else:
            non_nucleotide_chars_in_line = set(line) - VALID_NUCLEOTIDES
            if len(non_nucleotide_chars_in_line) > 0:
                msg = '{file}: Line {line} contains the following non-nucleotide characters: {chars}'.format(
                    file=filepath,
                    line=line_count,
                    chars=', '.join([str(x) for x in non_nucleotide_chars_in_line]))
                logging.warning(msg)
            seqs.append(line.upper())
        line_count += 1
    yield header, ''.join(seqs)


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
        with open(filepath, 'rU') as f:
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
