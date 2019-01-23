# -*- coding: utf-8 -*-

import logging
from typing import List, Any, Optional, Tuple, Union

import os
import re
from collections import defaultdict

from .const import SCHEME_FASTAS, REGEX_FASTQ, REGEX_FASTA
from .subtyping_params import SubtypingParams


def does_file_exist(filepath: str, force: bool):
    if filepath and os.path.exists(filepath):
        if not force:
            file_exists_err_fmt = 'File "{}" already exists! If you want to overwrite this output file run with opt "--force"'
            raise OSError(file_exists_err_fmt.format(filepath))
        else:
            logging.warning('File "{}" already exists, overwriting with "--force" - uh oh :S')


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
    for x, y in zip(a, b):
        if x != y:
            return False
    return True


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
    if len(fastqs) > 0:
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


NT_SUB = {x: y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}


def revcomp(s):
    """Reverse complement nucleotide sequence

    Args:
        s (str): nucleotide sequence

    Returns:
        str: reverse complement of `s` nucleotide sequence
    """
    return ''.join([NT_SUB[c] for c in s[::-1]])


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
        if args.max_missing_tiles:
            subtyping_params.max_perc_missing_tiles = args.max_missing_tiles
        if args.min_ambiguous_tiles:
            subtyping_params.min_ambiguous_tiles = args.min_ambiguous_tiles
        if args.max_intermediate_tiles:
            subtyping_params.max_perc_intermediate_tiles = args.max_intermediate_tiles
        if args.low_cov_warning:
            subtyping_params.min_coverage_warning = args.low_cov_warning
        if args.min_kmer_freq:
            subtyping_params.min_kmer_freq = args.min_kmer_freq
        if args.max_kmer_freq:
            subtyping_params.max_kmer_freq = args.max_kmer_freq

    return subtyping_params
