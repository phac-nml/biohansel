import logging
from typing import List, Tuple, Optional

import attr
import click
import pandas as pd

from biohansel.subtype import Subtype, subtype_contigs_samples, subtype_reads_samples
from biohansel.subtype.const import JSON_EXT_TMPL, SUBTYPE_SUMMARY_COLS
from biohansel.subtype.metadata import read_metadata_table, merge_metadata_with_summary_results
from biohansel.subtype.subtype_stats import subtype_counts
from biohansel.subtype.util import get_scheme_fasta, init_subtyping_params
from biohansel.utils import init_console_logger, does_file_exist, collect_inputs

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
@click.option('-v', '--verbose', count=True,
              help="Logging verbosity (-v for logging warnings; -vvv for logging debug info)")
def cli(verbose):
    """Find the closest matching NCBI RefSeq genomes or the genomes contained in your contigs or reads.
    """
    lvl = init_console_logger(verbose)
    logging.debug('Initialized logging with %s level', lvl)


@cli.command()
@click.option('-s', '--scheme', default='heidelberg',
              help='Scheme to use for subtyping (built-in: "heidelberg", "enteritidis"; '
                   'OR user-specified: /path/to/user/scheme)')
@click.option('--scheme-name', help='Custom user-specified SNP substyping scheme name')
@click.option('-M', '--scheme-metadata',
              help='Scheme subtype metadata table (CSV or tab-delimited format; must contain "subtype" column)')
@click.option('-p', '--paired-reads',
              # metavar=('forward_reads', 'reverse_reads'),
              type=(str, str),
              multiple=True,
              help='FASTQ paired-end reads')
@click.option('-i', '--input-fasta-genome-name',
              # metavar=('fasta_path', 'genome_name'),
              type=(str, str),
              multiple=True,
              help='fasta file path to genome name pair')
@click.option('-D', '--input-directory',
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='directory of input fasta files (.fasta|.fa|.fna) or FASTQ files (paired FASTQ should '
                   'have same basename with "_\d\.(fastq|fq)" postfix to be automatically paired) '
                   '(files can be Gzipped)')
@click.option('-o', '--output-summary',
              help='Subtyping summary output path (tab-delimited)')
@click.option('-O', '--output-tile-results',
              help='Subtyping tile matching output path (tab-delimited)')
@click.option('-S', '--output-simple-summary',
              help='Subtyping simple summary output path')
@click.option('--force',
              is_flag=True,
              help='Force existing output files to be overwritten')
@click.option('--json-output',
              is_flag=True,
              help='Output JSON representation of output files?')
@click.option('--min-kmer-freq',
              default=None,
              type=int,
              help='Min k-mer freq/coverage')
@click.option('--max-kmer-freq',
              default=None,
              type=int,
              help='Max k-mer freq/coverage')
@click.option('--low-coverage-threshold',
              default=None,
              type=float,
              help='Frequencies below this threshold are considered low coverage')
@click.option('--max-missing-tiles',
              default=None,
              type=float,
              help='Decimal proportion of maximum allowable missing tiles before being considered an error. (0.0 - 1.0)')
@click.option('--min-ambiguous-tiles',
              default=None,
              type=int,
              help='Minimum number of missing tiles to be considered an ambiguous result')
@click.option('--low-coverage-warning',
              default=None,
              type=int,
              help='Overall tile coverage below this value will trigger a low coverage warning')
@click.option('--max-intermediate-tiles',
              default=None,
              type=float,
              help='Decimal proportion of maximum allowable missing tiles to be considered an intermediate subtype. (0.0 - 1.0)')
@click.option('-t', '--threads',
              type=int,
              default=1,
              help='Number of parallel threads to run analysis (default=1)')
@click.argument('files', type=click.Path(exists=True), nargs=-1)
def subtype(scheme,
            scheme_name,
            scheme_metadata,
            paired_reads,
            input_fasta_genome_name,
            input_directory,
            output_summary,
            output_tile_results,
            output_simple_summary,
            force,
            json_output,
            min_kmer_freq,
            max_kmer_freq,
            low_coverage_threshold,
            max_missing_tiles,
            min_ambiguous_tiles,
            low_coverage_warning,
            max_intermediate_tiles,
            threads,
            files):
    output_summary_path = output_summary
    output_tile_results = output_tile_results
    output_simple_summary_path = output_simple_summary
    does_file_exist(output_simple_summary_path, force)
    does_file_exist(output_summary_path, force)
    does_file_exist(output_tile_results, force)
    scheme = scheme  # type: str
    scheme_name = scheme_name  # type: Optional[str]
    scheme_fasta = get_scheme_fasta(scheme)
    scheme_subtype_counts = subtype_counts(scheme_fasta)

    subtyping_params = init_subtyping_params(**locals())
    input_contigs, input_reads = collect_inputs(**locals())
    if len(input_contigs) == 0 and len(input_reads) == 0:
        raise Exception('No input files specified!')
    df_md = None
    if scheme_metadata:
        df_md = read_metadata_table(scheme_metadata)
    n_threads = threads

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
        if json_output:
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
            if json_output:
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
        if json_output:
            df_simple_summary.to_json(JSON_EXT_TMPL.format(output_simple_summary_path), **kwargs_for_pd_to_json)
