# -*- coding: utf-8 -*-
"""
Common utility functions and constants. 
"""
import logging
from typing import List, Any, Tuple, Union

import os
import re
from collections import defaultdict, Counter

LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

# Nucleotide substitution dictionary for reverse complementing a nucleotide sequence
NUCLEOTIDE_SUBSTITUTION = {x: y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}

REGEX_FASTQ = re.compile(r'^(.+)\.(fastq|fq|fastqsanger)(\.gz)?$')
REGEX_FASTA = re.compile(r'^.+\.(fasta|fa|fna|fas)(\.gz)?$')


def does_file_exist(filepath: str, force: bool) -> None:
    """Does a file exist? If so raise an OSError unless force is True

    Args:
        filepath: File path
        force: Overwrite file? Okay if the file already exists?

    Raises:
        OSError: if file already exists and it's not okay to have it be potentially overwritten.
    """
    if filepath and os.path.exists(filepath):
        if not force:
            raise OSError(
                f'File "{filepath}" already exists! If you want to overwrite this output file run with opt "--force"')
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


def _compare_subtypes(subtypes_a: List[Any], subtypes_b: List[Any]) -> bool:
    """Are 2 subtypes consistent with each other?

    Examples:
        ```python
        assert compare_subtypes([1,2,3], [1,2,3]) == True
        assert compare_subtypes([1,2,3], [1,2]) == True
        assert compare_subtypes([1], [1,2,3,4,5]) == True
        assert compare_subtypes([1], [2]) == False
        assert compare_subtypes([1,2,3,4], [1,2,3,5]) == False
        assert compare_subtypes([2,2], [1,2]) == False
        ```

    Args:
        subtypes_a: One list of subtype integers to compare
        subtypes_b:  Other list of subtype integers to compare

    Returns:
        True if lists are equal up to the last index of the smallest list in the comparison; False otherwise
    """
    for subtype_a, subtype_b in zip(subtypes_a, subtypes_b):
        if subtype_a != subtype_b:
            return False
    return True


def _compare_all_subtypes(subtypes: List[List[int]]) -> List[Tuple[List[int], List[int]]]:
    """Compare all subtypes against each other and return the pairs of inconsistent subtypes.

    Args:
        subtypes: Subtypes to compare against each other

    Returns:
        Pairs of subtypes that are inconsistent with each other.
    """
    inconsistent_subtypes = []
    for i in range(len(subtypes) - 1):
        subtypes_a = subtypes[i]
        for j in range(i + 1, len(subtypes)):
            subtypes_b = subtypes[j]
            is_consistent = _compare_subtypes(subtypes_a, subtypes_b)
            if not is_consistent:
                inconsistent_subtypes.append((subtypes_a, subtypes_b))
    return inconsistent_subtypes


def _inconsistent_subtype_pairs_to_list(inconsistent_subtype_pairs: List[Tuple[List[int], List[int]]]) -> List[str]:
    potentially_inconsistent_subtypes = []
    for subtypes_a, subtypes_b in inconsistent_subtype_pairs:
        subtype_a = '.'.join([str(x) for x in subtypes_a])
        subtype_b = '.'.join([str(x) for x in subtypes_b])
        potentially_inconsistent_subtypes += [subtype_a, subtype_b]
    return potentially_inconsistent_subtypes


def find_inconsistent_subtypes(subtypes: List[List[int]]) -> List[str]:
    """Find all inconsistent subtypes within a given list of subtypes.

    Args:
        subtypes: List of subtypes to compare against each other

    Returns:
        List of inconsistent subtype strings
    """
    inconsistent_subtype_pairs = _compare_all_subtypes(subtypes)
    potentially_inconsistent_subtypes = _inconsistent_subtype_pairs_to_list(inconsistent_subtype_pairs)
    return [subtype for subtype, freq in Counter(potentially_inconsistent_subtypes).most_common() if freq >= 1]


def collect_fastq_from_dir(input_directory: str) -> List[Union[str, Tuple[List[str], str]]]:
    fastqs = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and REGEX_FASTQ.match(x):
            fastqs.append(full_file_path)
    if len(fastqs) > 0:
        logging.info(f'Found {len(fastqs)} FASTQ files in "{input_directory}"')
        reads_from_dir = group_fastqs(fastqs)
        logging.info(f'Collected {len(reads_from_dir)} read sets from {len(fastqs)} FASTQ files in "{input_directory}"')
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


def revcomp(nucleotide_sequence: str) -> str:
    """Reverse complement nucleotide sequence

    Args:
        nucleotide_sequence: nucleotide sequence

    Returns:
        Reverse complemented nucleotide sequence
    """
    return ''.join([NUCLEOTIDE_SUBSTITUTION[nt] for nt in nucleotide_sequence[::-1]])


def is_gzipped(filepath: str) -> bool:
    """Is a file gzipped? Does the file path end with `.gz`?"""
    return bool(re.match(r'^.+\.gz$', filepath))


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)


def collect_inputs(files: List[str],
                   input_fasta_genome_name: List[Tuple[str, str]],
                   input_directory: str,
                   paired_reads: List[Union[List[str], Tuple[str, str]]],
                   **kwargs) -> Tuple[List[Tuple[str, str]], List[Tuple[List[str], str]]]:
    """Collect all input files for analysis

    Sample names are derived from the base filename with no extensions.
    Sequencing reads are paired if they share a common filename name without "_\d".
    Filepaths for contigs and reads files are collected from an input directory if provided.

    Args:
        args: ArgumentParser.parse_args() output

    Returns:
        List of (contig filename, sample name)
        List of ([reads filepaths], sample name)
    """
    input_genomes = []
    reads = []
    if files:
        fastas = [x for x in files if REGEX_FASTA.match(x)]
        fastqs = [x for x in files if REGEX_FASTQ.match(x)]
        if len(fastas) > 0:
            for fasta_path in fastas:
                fasta_path = os.path.abspath(fasta_path)
                if os.path.exists(fasta_path):
                    genome_name = genome_name_from_fasta_path(fasta_path)
                    input_genomes.append((fasta_path, genome_name))
                else:
                    logging.error(f'Input fasta "{fasta_path}" does not exist!')
        if len(fastqs) > 0:
            grouped_fastqs = group_fastqs(fastqs)
            logging.info(f'Grouped {len(fastqs)} fastqs into {len(grouped_fastqs)} groups')
            reads += grouped_fastqs
    if input_fasta_genome_name:
        for fasta_path, genome_name in input_fasta_genome_name:
            input_genomes.append((os.path.abspath(fasta_path), genome_name))
    if input_directory:
        logging.info(f'Searching dir "{input_directory}" for FASTA files')
        input_genomes += collect_fasta_from_dir(input_directory)
        logging.info(f'Searching dir "{input_directory}" for FASTQ files')
        reads += collect_fastq_from_dir(input_directory)
    if paired_reads:
        for filepaths in paired_reads:
            if not isinstance(filepaths, (list, tuple)):
                logging.warning(f'Paired end reads not list or tuple {filepaths}')
                continue
            filenames = [os.path.basename(filepath) for filepath in filepaths]
            common_prefix = os.path.commonprefix(filenames)
            genome_name = re.sub(r'[\W\_]+$', r'', common_prefix)
            if genome_name == '':
                genome_name = filenames[0]
            reads.append((filepaths, genome_name))
    return input_genomes, reads
