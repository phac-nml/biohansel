# -*- coding: utf-8 -*-

import logging
import os
import re
from collections import defaultdict
from itertools import product
from typing import List, Any, Optional, Tuple, Union

import pandas as pd

from .const import SCHEME_FASTAS, REGEX_FASTQ, REGEX_FASTA, bases_dict
from .parsers import parse_fasta
from .subtyping_params import SubtypingParams


def does_file_exist(filepath: str, force: bool):
    if filepath and os.path.exists(filepath):
        if not force:
            msg = f'File "{filepath}" already exists! If you want to overwrite this output file run with opt "--force"'
            raise OSError(msg)
        else:
            logging.warning(f'File "{filepath}" already exists, overwriting with "--force"')


def genome_name_from_fasta_path(fasta_path: str) -> str:
    """Extract genome name from fasta filename

    Get the filename without directory and remove the file extension.

    Example:
        With fasta file path ``/path/to/genome_1.fasta``::

            fasta_path = '/path/to/genome_1.fasta'
            genome_name = genome_name_from_fasta_path(fasta_path)
            print(genome_name)
            # => "genome_1"

    Args:
        fasta_path (str): fasta file path

    Returns:
        str: genome name
    """
    filename = os.path.basename(fasta_path)
    filename = re.sub(r'\.gz$', '', filename)
    return re.sub(r'\.(fa|fas|fasta|fna|\w{1,})(\.gz)?$', '', filename)


def compare_subtypes(a: List[Any], b: List[Any]) -> bool:
    return all(x == y for x, y in zip(a, b))


def find_inconsistent_subtypes(subtypes: List[List[int]]) -> List[str]:
    from collections import Counter
    incon = []
    for i in range(len(subtypes) - 1):
        a = subtypes[i]
        for j in range(i + 1, len(subtypes)):
            b = subtypes[j]
            is_consistent = compare_subtypes(a, b)
            if not is_consistent:
                incon.append((a, b))
    l = []
    for a, b in incon:
        astr = '.'.join([str(x) for x in a])
        bstr = '.'.join([str(x) for x in b])
        l += [astr, bstr]
    c = Counter(l)
    incon_subtypes = []
    for subtype, freq in c.most_common():
        if freq >= 1:
            incon_subtypes.append(subtype)
        else:
            break
    return incon_subtypes


def get_scheme_fasta(scheme: str) -> str:
    if scheme in SCHEME_FASTAS:
        scheme_fasta = SCHEME_FASTAS[scheme]['file']
    elif os.path.exists(scheme) and os.path.isfile(scheme):
        scheme_fasta = scheme
    else:
        raise FileNotFoundError('Could not find user-specified subtyping scheme fasta "%s"', scheme)
    return scheme_fasta


def get_scheme_params(scheme: str) -> Optional[SubtypingParams]:
    if scheme in SCHEME_FASTAS:
        return SCHEME_FASTAS[scheme]['subtyping_params']


def get_scheme_version(scheme: str) -> Optional[str]:
    if scheme in SCHEME_FASTAS:
        version = SCHEME_FASTAS[scheme]['version']  # type: str
        return version
    return None


def collect_fastq_from_dir(input_directory: str) -> List[Union[str, Tuple[List[str], str]]]:
    fastqs = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and REGEX_FASTQ.match(x):
            fastqs.append(full_file_path)
    if fastqs:
        logging.info('Found %s FASTQ files in %s',
                     len(fastqs),
                     input_directory)
        reads_from_dir = group_fastqs(fastqs)
        logging.info('Collected %s read sets from %s FASTQ files in %s',
                     len(reads_from_dir),
                     len(fastqs),
                     input_directory)
        return reads_from_dir
    return []


def group_fastqs(fastqs: List[str]) -> List[Tuple[List[str], str]]:
    """Group FASTQs based on common base filename

    For example, if you have 2 FASTQs:

    - reads_1.fastq
    - reads_2.fastq

    The common name would be `reads` and the files would be grouped based on that common name.

    Args:
        fastqs: FASTQ file paths

    Returns:
        list of grouped FASTQs grouped by common base filename
    """
    genome_fastqs = defaultdict(list)
    for fastq in fastqs:
        filename = os.path.basename(fastq)
        basefilename = re.sub(r'_\d', '', REGEX_FASTQ.sub(r'\1', filename))
        genome_fastqs[basefilename].append(fastq)
    return [(fastq_paths, genome_name) for genome_name, fastq_paths in genome_fastqs.items()]


def collect_fasta_from_dir(input_directory: str) -> List[Tuple[str, str]]:
    input_genomes = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and REGEX_FASTA.match(x):
            genome_name = genome_name_from_fasta_path(full_file_path)
            input_genomes.append((full_file_path, genome_name))
    return input_genomes


NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')


def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]


def is_gzipped(p: str) -> bool:
    return bool(re.match(r'^.+\.gz$', p))


def init_subtyping_params(args: Optional[Any] = None,
                          scheme: Optional[str] = None) -> SubtypingParams:
    """Initialize subtyping parameters based on command-line arguments and scheme defaults

    Args:
        args: ArgumentParser.parse_args() output
        scheme: Scheme name e.g. "heidelberg"

    Returns:
        SubtypingParams with user-supplied values then scheme defaults then global defaults loaded
    """
    subtyping_params = get_scheme_params(scheme)
    if subtyping_params is None:
        subtyping_params = SubtypingParams()
    if args is not None:
        if args.low_cov_depth_freq:
            subtyping_params.low_coverage_depth_freq = args.low_cov_depth_freq
        if args.max_missing_kmers:
            subtyping_params.max_perc_missing_kmers = args.max_missing_kmers
        if args.min_ambiguous_kmers:
            subtyping_params.min_ambiguous_kmers = args.min_ambiguous_kmers
        if args.max_intermediate_kmers:
            subtyping_params.max_perc_intermediate_kmers = args.max_intermediate_kmers
        if args.low_cov_warning:
            subtyping_params.min_coverage_warning = args.low_cov_warning
        if args.min_kmer_freq:
            subtyping_params.min_kmer_freq = args.min_kmer_freq
        if args.min_kmer_frac:
            subtyping_params.min_kmer_frac = args.min_kmer_frac
        if args.max_kmer_freq:
            subtyping_params.max_kmer_freq = args.max_kmer_freq
        if args.max_degenerate_kmers:
            subtyping_params.max_degenerate_kmers = args.max_degenerate_kmers

    return subtyping_params


def df_field_fillna(df: pd.DataFrame, field: str = 'subtype', na: str = '#N/A') -> pd.DataFrame:
    df[field].replace('', na, inplace=True)
    df[field].fillna(value=na, inplace=True)
    df[field] = df[field].astype(str)
    return df


def check_total_kmers(scheme_fasta, max_degenerate_kmers):
    """Checks that the number of kmers about to be created is not at too high a computation or time cost

    Args:
         scheme_fasta: Kmer sequences from the SNV scheme
         max_degenerate_kmers:  The max kmers allowed by the scheme

    Raises:
        ValueError if number of created kmers is greater than the max degenerate kmers argument
    """
    kmer_number = 0
    for header, sequence in parse_fasta(scheme_fasta):
        value = 1
        for char in sequence:
            length_key = len(bases_dict[char])
            value *= length_key
        kmer_number += value
    if kmer_number * 2 > max_degenerate_kmers:
        raise ValueError(f'Your current scheme contains "{kmer_number * 2}" '
                         f'kmers which is over the current max degenerate '
                         f'kmers check of "{max_degenerate_kmers}" '
                         f'(Maximum recommended k-mers is 100,000). '
                         f'It is not advised to run this scheme due to the '
                         f'time and memory usage required to give an output '
                         f'with this many kmers loaded. If you still want to '
                         f'run this scheme, add the command line check of '
                         f'"--max-degenerate-kmers {kmer_number * 2 + 1}" '
                         f'at the end of your previous command.')


def expand_degenerate_bases(seq):
    """List all possible kmers for a scheme given a degenerate base

    Args:
         Scheme_kmers from SNV scheme fasta file

    Returns:
         List of all possible kmers given a degenerate base or not
    """

    return list(map("".join, product(*map(bases_dict.get, seq))))