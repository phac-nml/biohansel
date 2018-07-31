#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys
import re
import os
from typing import Optional, List, Any, Tuple

import attr
import pandas as pd

from . import program_desc, __version__
from .const import SUBTYPE_SUMMARY_COLS, REGEX_FASTQ, REGEX_FASTA, JSON_EXT_TMPL
from .subtype import Subtype
from .subtype_stats import subtype_counts
from .subtyper import \
    subtype_contigs_samples, \
    subtype_reads_samples
from .metadata import read_metadata_table, merge_metadata_with_summary_results
from .utils import \
    genome_name_from_fasta_path, \
    get_scheme_fasta, \
    does_file_exist, \
    collect_fastq_from_dir, \
    group_fastqs, \
    collect_fasta_from_dir, \
    init_subtyping_params

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
                        help='Input genome FASTA/FASTQ files (can be Gzipped)')
    parser.add_argument('-s', '--scheme',
                        default='heidelberg',
                        help='Scheme to use for subtyping (built-in: "heidelberg", "enteritidis"; OR user-specified: '
                             '/path/to/user/scheme)')
    parser.add_argument('--scheme-name',
                        help='Custom user-specified SNP substyping scheme name')
    parser.add_argument('-M', '--scheme-metadata',
                        help='Scheme subtype metadata table (CSV or tab-delimited format; must contain "subtype" column)')
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
                        help='directory of input fasta files (.fasta|.fa|.fna) or FASTQ files (paired FASTQ should '
                             'have same basename with "_\d\.(fastq|fq)" postfix to be automatically paired) '
                             '(files can be Gzipped)')
    parser.add_argument('-o', '--output-summary',
                        help='Subtyping summary output path (tab-delimited)')
    parser.add_argument('-O', '--output-tile-results',
                        help='Subtyping tile matching output path (tab-delimited)')
    parser.add_argument('-S', '--output-simple-summary',
                        help='Subtyping simple summary output path')
    parser.add_argument('--force',
                        action='store_true',
                        help='Force existing output files to be overwritten')
    parser.add_argument('--json',
                        action='store_true',
                        help='Output JSON representation of output files')
    parser.add_argument('--min-kmer-freq',
                        type=int,
                        help='Min k-mer freq/coverage')
    parser.add_argument('--max-kmer-freq',
                        type=int,
                        help='Max k-mer freq/coverage')
    parser.add_argument('--low-cov-depth-freq',
                        type=int,
                        help='Frequencies below this coverage are considered low coverage')
    parser.add_argument('--max-missing-tiles',
                        type=float,
                        help='Decimal proportion of maximum allowable missing tiles before being considered an error. (0.0 - 1.0)')
    parser.add_argument('--min-ambiguous-tiles',
                        type=int,
                        help='Minimum number of missing tiles to be considered an ambiguous result')
    parser.add_argument('--low-cov-warning',
                        type=int,
                        help='Overall tile coverage below this value will trigger a low coverage warning')
    parser.add_argument('--max-intermediate-tiles',
                        type=float,
                        help='Decimal proportion of maximum allowable missing tiles to be considered an intermediate subtype. (0.0 - 1.0)')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run analysis (default=1)')
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s {}'.format(__version__))
    return parser


def collect_inputs(args: Any) -> Tuple[List[Tuple[str, str]], List[Tuple[List[str], str]]]:
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
    if args.files:
        fastas = [x for x in args.files if REGEX_FASTA.match(x)]
        fastqs = [x for x in args.files if REGEX_FASTQ.match(x)]
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
    return input_genomes, reads


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
    does_file_exist(output_simple_summary_path, args.force)
    does_file_exist(output_summary_path, args.force)
    does_file_exist(output_tile_results, args.force)
    scheme = args.scheme  # type: str
    scheme_name = args.scheme_name  # type: Optional[str]
    scheme_fasta = get_scheme_fasta(scheme)
    scheme_subtype_counts = subtype_counts(scheme_fasta)
    logging.debug(args)
    subtyping_params = init_subtyping_params(args, scheme)
    input_contigs, input_reads = collect_inputs(args)
    if len(input_contigs) == 0 and len(input_reads) == 0:
        raise Exception('No input files specified!')
    df_md = None
    if args.scheme_metadata:
        df_md = read_metadata_table(args.scheme_metadata)
    n_threads = args.threads

    subtype_results = []  # type: List[Tuple[Subtype, pd.DataFrame]]
    if len(input_contigs) > 0:
        contigs_results = subtype_contigs_samples(input_genomes=input_contigs,
                                                  scheme=scheme,
                                                  scheme_name=scheme_name,
                                                  subtyping_params=subtyping_params,
                                                  scheme_subtype_counts=scheme_subtype_counts,
                                                  n_threads=n_threads)
        logging.info('Generated %s subtyping results from %s contigs samples', len(contigs_results), len(input_contigs))
        subtype_results += contigs_results
    if len(input_reads) > 0:
        reads_results = subtype_reads_samples(reads=input_reads,
                                              scheme=scheme,
                                              scheme_name=scheme_name,
                                              subtyping_params=subtyping_params,
                                              scheme_subtype_counts=scheme_subtype_counts,
                                              n_threads=n_threads)
        logging.info('Generated %s subtyping results from %s contigs samples', len(reads_results), len(input_reads))
        subtype_results += reads_results

    dfs = [df for st, df in subtype_results]  # type: List[pd.DataFrame]
    dfsummary = pd.DataFrame([attr.asdict(st) for st, df in subtype_results])
    dfsummary = dfsummary[SUBTYPE_SUMMARY_COLS]

    if dfsummary['avg_tile_coverage'].isnull().all():
        dfsummary = dfsummary.drop(labels='avg_tile_coverage', axis=1)

    if df_md is not None:
        dfsummary = merge_metadata_with_summary_results(dfsummary, df_md)

    kwargs_for_pd_to_table = dict(sep='\t', index=None, float_format='%.3f')
    kwargs_for_pd_to_json = dict(orient='records')

    if output_summary_path:
        dfsummary.to_csv(output_summary_path, **kwargs_for_pd_to_table)
        if args.json:
            dfsummary.to_json(JSON_EXT_TMPL.format(output_summary_path), **kwargs_for_pd_to_json)
        logging.info('Wrote subtyping output summary to %s', output_summary_path)
    else:
        # if no output path specified for the summary results, then print to stdout
        print(dfsummary.to_csv(sep='\t', index=None))

    if output_tile_results:
        if len(dfs) > 0:
            dfall = pd.concat(dfs)  # type: pd.DataFrame
            dfall.to_csv(output_tile_results, **kwargs_for_pd_to_table)
            logging.info('Tile results written to "{}".'.format(output_tile_results))
            if args.json:
                dfall.to_json(JSON_EXT_TMPL.format(output_tile_results), **kwargs_for_pd_to_json)
                logging.info(
                    'Tile results written to "{}" in JSON format.'.format(JSON_EXT_TMPL.format(output_tile_results)))
        else:
            logging.error(
                'No tile results generated. No tile results file written to "{}".'.format(output_tile_results))

    if output_simple_summary_path:
        if 'avg_tile_coverage' in dfsummary.columns:
            df_simple_summary = dfsummary[['sample', 'subtype', 'avg_tile_coverage', 'qc_status', 'qc_message']]
        else:
            df_simple_summary = dfsummary[['sample', 'subtype', 'qc_status', 'qc_message']]

        if df_md is not None:
            df_simple_summary = merge_metadata_with_summary_results(df_simple_summary, df_md)

        df_simple_summary.to_csv(output_simple_summary_path, **kwargs_for_pd_to_table)
        if args.json:
            df_simple_summary.to_json(JSON_EXT_TMPL.format(output_simple_summary_path), **kwargs_for_pd_to_json)


if __name__ == '__main__':
    main()
