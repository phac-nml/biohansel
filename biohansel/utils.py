# -*- coding: utf-8 -*-

import logging
from typing import List, Any, Tuple, Union

import os
import re
from collections import defaultdict

NT_SUB = {x: y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}
REGEX_FASTQ = re.compile(r'^(.+)\.(fastq|fq|fastqsanger)(\.gz)?$')
REGEX_FASTA = re.compile(r'^.+\.(fasta|fa|fna|fas)(\.gz)?$')


def does_file_exist(filepath: str, force: bool):
    if filepath and os.path.exists(filepath):
        if not force:
            raise OSError(f'File "{filepath}" already exists! If you want to overwrite this output file run with opt "--force"')
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


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)


def collect_inputs(files,
                   input_fasta_genome_name,
                   input_directory,
                   paired_reads,
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
            logging.info('# of input fastas %s', len(fastas))
            for fasta_path in fastas:
                fasta_path = os.path.abspath(fasta_path)
                if os.path.exists(fasta_path):
                    genome_name = genome_name_from_fasta_path(fasta_path)
                    input_genomes.append((fasta_path, genome_name))
                else:
                    logging.error('Input fasta "%s" does not exist!', fasta_path)
        if len(fastqs) > 0:
            logging.info('# of input fastqs %s', len(fastqs))
            grouped_fastqs = group_fastqs(fastqs)
            logging.info('Grouped %s fastqs into %s groups',
                         len(fastqs),
                         len(grouped_fastqs))
            reads += grouped_fastqs
    if input_fasta_genome_name:
        for fasta_path, genome_name in input_fasta_genome_name:
            input_genomes.append((os.path.abspath(fasta_path), genome_name))
    if input_directory:
        logging.info('Searching dir "%s" for FASTA files', input_directory)
        input_genomes += collect_fasta_from_dir(input_directory)
        logging.info('Searching dir "%s" for FASTQ files', input_directory)
        reads += collect_fastq_from_dir(input_directory)
    if paired_reads:
        for x in paired_reads:
            if not isinstance(x, (list, tuple)):
                logging.warning('Paired end reads not list or tuple %s', x)
                continue
            filenames = [os.path.basename(y) for y in x]
            common_prefix = os.path.commonprefix(filenames)
            genome_name = re.sub(r'[\W\_]+$', r'', common_prefix)
            if genome_name == '':
                genome_name = filenames[0]
            reads.append((x, genome_name))
    return input_genomes, reads


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
