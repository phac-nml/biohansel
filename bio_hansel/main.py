#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from typing import Optional

import attr
import logging
import re
import pandas as pd
from collections import defaultdict

from bio_hansel import program_name, program_desc, __version__
from bio_hansel.const import SUBTYPE_SUMMARY_COLS
from bio_hansel.subtyper import subtype_fasta, subtype_reads
from bio_hansel.subtype_stats import subtype_counts
from bio_hansel.subtyping_params import SubtypingParams
from bio_hansel.utils import genome_name_from_fasta_path, get_scheme_fasta, out_files_exists, get_scheme_params

SCRIPT_NAME = 'hansel'
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'


def init_console_logger(logging_verbosity=3):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]

    logging.basicConfig(format=LOG_FORMAT, level=lvl)


def init_parser():
    parser = argparse.ArgumentParser(prog=SCRIPT_NAME,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=program_desc)
    parser.add_argument('files',
                        metavar='F',
                        nargs='*',
                        help='Input genome FASTA/FASTQ files')
    parser.add_argument('-s', '--scheme',
                        default='heidelberg',
                        help='Scheme to use for subtyping (built-in: "heidelberg", "enteritidis"; OR user-specified: /path/to/user/scheme)')
    parser.add_argument('--scheme-name',
                        help='Custom user-specified SNP substyping scheme name')
    parser.add_argument('-p', '--paired-reads',
                        nargs=2,
                        metavar=('forward_reads', 'reverse_reads'),
                        action='append',
                        help='FASTQ paired-end reads')
    parser.add_argument('-i', '--input-fasta-genome-name',
                        nargs=2,
                        metavar=('fasta_path', 'genome_name'),
                        action='append',
                        help='fasta file path to genome name pair')
    parser.add_argument('-D', '--input-directory',
                        help='directory of input fasta files (.fasta|.fa|.fna) or FASTQ files (paired FASTQ should have same basename with "_\d\.(fastq|fq)" postfix to be automatically paired)')
    parser.add_argument('-o', '--output-summary',
                        help='Subtyping summary output path (tab-delimited)')
    parser.add_argument('-O', '--output-tile-results',
                        help='Subtyping tile matching output path (tab-delimited)')
    parser.add_argument('-S', '--output-simple-summary',
                        help='Subtyping simple summary output path')
    parser.add_argument('--force',
                        action='store_true',
                        help='Force existing output files to be overwritten')
    parser.add_argument('--min-kmer-freq',
                        type=int,
                        help='Min k-mer freq/coverage')
    parser.add_argument('--max-kmer-freq',
                        type=int,
                        help='Max k-mer freq/coverage')
    # Changes
    parser.add_argument('--low-cov-depth-freq',
                        type=int,
                        help='Frequencies below this coverage are considered low coverage')
    parser.add_argument('--max-missing-tiles',
                        type=float,
                        help='Decimal proportion of maximum allowable missing tiles before being considered an error. (0.0 - 1.0)')
    parser.add_argument('--min-ambiguous-tiles',
                        type=int,
                        help='Minimum number of missing tiles to be considered an ambiguous result')
    parser.add_argument('--max-intermediate-tiles',
                        type=float,
                        help='Decimal proportion of maximum allowable missing tiles to be considered an intermediate subtype. (0.0 - 1.0)')
    # Changes
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run analysis (default=1)')
    parser.add_argument('-T', '--tmp-dir',
                        default='/tmp',
                        help='Base temporary working directory for intermediate analysis files')
    parser.add_argument('-K', '--keep-tmp',
                        default=False,
                        action='store_true',
                        help='Keep temporary analysis files')
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


def main():
    parser = init_parser()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    output_summary_path = args.output_summary
    output_tile_results = args.output_tile_results
    output_simple_summary_path = args.output_simple_summary
    out_files_exists(output_simple_summary_path, args.force)
    out_files_exists(output_summary_path, args.force)
    out_files_exists(output_tile_results, args.force)
    scheme = args.scheme  # type: str
    scheme_name = args.scheme_name  # type: Optional[str]
    scheme_fasta = get_scheme_fasta(scheme)
    scheme_subtype_counts = subtype_counts(scheme_fasta)
    input_genomes = []
    reads = []
    logging.debug(args)

    subtyping_params = get_scheme_params(scheme)
    if not subtyping_params:
        subtyping_params = SubtypingParams()

    if args.low_cov_depth_freq:
        subtyping_params.low_coverage_depth_freq = args.low_cov_depth_freq
    if args.max_missing_tiles:
        subtyping_params.max_perc_missing_tiles = args.max_missing_tiles
    if args.min_ambiguous_tiles:
        subtyping_params.min_ambiguous_tiles = args.min_ambiguous_tiles
    if args.max_intermediate_tiles:
        subtyping_params.max_perc_intermediate_tiles = args.max_intermediate_tiles

    if args.files:
        fastas = [x for x in args.files if re.match(r'^.+\.(fasta|fa|fna)$', x)]
        fastqs = [x for x in args.files if re.match(r'^.+\.(fastq|fq)$', x)]
        if fastas:
            logging.info('# of input fastas %s', len(fastas))
            for fasta_path in fastas:
                fasta_path = os.path.abspath(fasta_path)
                if os.path.exists(fasta_path):
                    genome_name = genome_name_from_fasta_path(fasta_path)
                    input_genomes.append((fasta_path, genome_name))
                else:
                    logging.error('Input fasta "%s" does not exist!', fasta_path)
        if fastqs:
            logging.info('# of input fastqs %s', len(fastqs))
            grouped_fastqs = group_fastqs(fastqs)
            logging.info('Grouped %s fastqs into %s groups',
                         len(fastqs),
                         len(grouped_fastqs))
            reads += grouped_fastqs

    if args.input_fasta_genome_name:
        for fasta_path, genome_name in args.input_fasta_genome_name:
            input_genomes.append((os.path.abspath(fasta_path), genome_name))

    if args.input_directory:
        logging.info('Searching dir "%s" for FASTA files', args.input_directory)
        input_genomes += collect_fasta_from_dir(args.input_directory)
        logging.info('Searching dir "%s" for FASTQ files', args.input_directory)
        reads += collect_fastq_from_dir(args.input_directory)

    if args.paired_reads:
        for x in args.paired_reads:
            if not isinstance(x, (list, tuple)):
                logging.warning('Paired end reads not list or tuple %s', x)
                continue
            filenames = [os.path.basename(y) for y in x]
            common_prefix = os.path.commonprefix(filenames)
            genome_name = re.sub(r'[\W\_]+$', r'', common_prefix)
            if genome_name == '':
                genome_name = filenames[0]
            reads.append((x, genome_name))

    if len(input_genomes) == 0 and len(reads) == 0:
        raise Exception('No input files specified!')

    n_threads = args.threads
    tmp_dir = args.tmp_dir

    subtype_results = []
    dfs = []
    if input_genomes:
        if n_threads == 1:
            logging.info('Serial single threaded run mode on %s input genomes', len(input_genomes))
            outputs = [subtype_fasta(subtyping_params,
                                     scheme,
                                     input_fasta,
                                     genome_name,
                                     tmp_dir=tmp_dir,
                                     scheme_name=scheme_name,
                                     scheme_subtype_counts=scheme_subtype_counts)
                       for input_fasta, genome_name in input_genomes]
        else:
            from multiprocessing import Pool
            logging.info('Initializing thread pool with %s threads', n_threads)
            pool = Pool(processes=n_threads)
            logging.info('Running analysis asynchronously on %s input genomes', len(input_genomes))
            res = [pool.apply_async(subtype_fasta, (subtyping_params,
                                                    scheme,
                                                    input_fasta,
                                                    genome_name,
                                                    tmp_dir,
                                                    scheme_name,
                                                    scheme_subtype_counts))
                   for input_fasta, genome_name in input_genomes]

            logging.info('Parallel analysis complete! Retrieving analysis results')
            outputs = [x.get() for x in res]

        for subtype, df in outputs:
            if df is not None:
                dfs.append(df)
            subtype_results.append(attr.asdict(subtype))

    if reads:
        outputs = [subtype_reads(subtyping_params,
                                 scheme=scheme,
                                 reads=r,
                                 genome_name=genome_name,
                                 tmp_dir=tmp_dir,
                                 threads=n_threads,
                                 scheme_name=scheme_name,
                                 scheme_subtype_counts=scheme_subtype_counts)
                   for r, genome_name in reads]

        for subtype, df in outputs:
            if df is not None:
                dfs.append(df)
            subtype_results.append(attr.asdict(subtype))

    dfall = pd.concat(dfs)  # type: pd.DataFrame
    dfsummary = pd.DataFrame(subtype_results)
    dfsummary = dfsummary[SUBTYPE_SUMMARY_COLS]

    df_simple_summary = dfsummary[['sample', 'subtype', 'qc_status', 'qc_message']]

    if output_summary_path:
        dfsummary.to_csv(output_summary_path, sep='\t', index=None)
        logging.info('Wrote subtyping output summary to %s', output_summary_path)
    else:
        print(dfsummary.to_csv(sep='\t', index=None))

    if output_tile_results:
        dfall.to_csv(output_tile_results, sep='\t', index=None)

    if output_simple_summary_path:
        df_simple_summary.to_csv(output_simple_summary_path, sep='\t', index=None)


def collect_fasta_from_dir(input_directory):
    input_genomes = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and re.match(r'^.+\.(fasta|fa|fna)$', x):
            genome_name = genome_name_from_fasta_path(full_file_path)
            input_genomes.append((full_file_path, genome_name))
    return input_genomes


def group_fastqs(fastqs):
    reads = []
    genome_fastqs = defaultdict(list)
    for fastq in fastqs:
        filename = os.path.basename(fastq)
        basefilename = re.sub(r'_\d\.(fastq|fq)$', '', filename)
        genome_fastqs[basefilename].append(fastq)
    for genome_name, fastq_paths in genome_fastqs.items():
        reads.append((fastq_paths, genome_name))
    return reads


def collect_fastq_from_dir(input_directory):
    fastqs = []
    for x in os.listdir(input_directory):
        full_file_path = os.path.abspath(os.path.join(input_directory, x))
        if os.path.isfile(full_file_path) and re.match(r'^.+\.(fastq|fq)$', x):
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


if __name__ == '__main__':
    main()
