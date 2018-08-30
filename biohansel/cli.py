import logging
import os
from typing import Union, Optional, Tuple, List

import attr
import click
import pandas as pd

from biohansel.create.display_tree import display_tree
from biohansel.create.cluster_generator import find_clusters
from biohansel.create.io.output import write_sequence_file
from biohansel.create.io.parsers import parse_vcf, parse_sequence_file
from biohansel.create.schema_generator import get_sequences, group_snvs
from biohansel.subtype import subtype_contigs_samples, subtype_reads_samples, Subtype
from biohansel.subtype.const import SUBTYPE_SUMMARY_COLS, JSON_EXT_TMPL
from biohansel.subtype.metadata import read_metadata_table, merge_metadata_with_summary_results
from biohansel.subtype.subtype_stats import subtype_counts
from biohansel.subtype.util import get_scheme_fasta, init_subtyping_params
from biohansel.utils import does_file_exist, collect_inputs, genome_name_from_fasta_path, init_console_logger
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
@click.option('-v', '--verbose', count=True,
              help="Logging verbosity (-v for logging warnings; -vvv for logging debug info)")
def cli(verbose):
    """Subtype with a biohansel scheme or create a scheme for your organism of interest
    """
    init_console_logger(verbose)
    logging.info(f'Initialized logging at "{logging.getLevelName(logging.getLogger().level)}" level')


def check_between_0_and_1_inclusive(ctx: click.Context,
                                    param: click.Option,
                                    value: Optional[Union[int, float]]) -> Optional[Union[int, float]]:
    if value is None or 0.0 <= value <= 1.0:
        return value
    else:
        raise click.BadParameter('value needs to be between 0.0 and 1.0 inclusive!')


def check_positive_number(ctx: click.Context,
                          param: click.Option,
                          value: Optional[Union[int, float]]) -> Optional[Union[int, float]]:
    if value is None or value >= 0:
        return value
    else:
        raise click.BadParameter('value must be 0 or greater!')


@cli.command()
@click.option('-s', '--scheme', default='heidelberg',
              help='Scheme to use for subtyping (built-in: "heidelberg", "enteritidis"; '
                   'OR user-specified: /path/to/user/scheme)')
@click.option('--scheme-name', help='Custom user-specified SNP substyping scheme name')
@click.option('-M', '--scheme-metadata',
              help='Scheme subtype metadata table (CSV or tab-delimited format; must contain "subtype" column)')
@click.option('-p', '--paired-reads',
              metavar='FORWARD_READS REVERSE_READS',
              type=(str, str),
              multiple=True,
              help='FASTQ paired-end reads')
@click.option('-i', '--input-fasta-genome-name',
              metavar='FASTA_PATH GENOME_NAME',
              type=(str, str),
              multiple=True,
              help='fasta file path to genome name pair')
@click.option('-D', '--input-directory',
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='directory of input fasta files (.fasta|.fa|.fna) or FASTQ files (paired FASTQ should '
                   'have same basename with "_\d\.(fastq|fq)" postfix to be automatically paired) '
                   '(files can be Gzipped)')
@click.option('-o', '--output-summary-path',
              help='Subtyping summary output path (tab-delimited)')
@click.option('-O', '--output-tile-results',
              help='Subtyping tile matching output path (tab-delimited)')
@click.option('-S', '--output-simple-summary-path',
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
              callback=check_positive_number,
              help='Min k-mer freq/coverage')
@click.option('--max-kmer-freq',
              default=None,
              type=int,
              callback=check_positive_number,
              help='Max k-mer freq/coverage')
@click.option('--low-coverage-threshold',
              default=None,
              type=float,
              callback=check_positive_number,
              help='Frequencies below this threshold are considered low coverage')
@click.option('--max-missing-tiles',
              default=None,
              type=float,
              callback=check_between_0_and_1_inclusive,
              help='Decimal proportion of maximum allowable missing tiles before being considered an error. (0.0 - 1.0)')
@click.option('--min-ambiguous-tiles',
              default=None,
              type=int,
              callback=check_positive_number,
              help='Minimum number of missing tiles to be considered an ambiguous result')
@click.option('--low-coverage-warning',
              default=None,
              type=int,
              callback=check_positive_number,
              help='Overall tile coverage below this value will trigger a low coverage warning')
@click.option('--max-intermediate-tiles',
              default=None,
              type=float,
              callback=check_between_0_and_1_inclusive,
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
            output_summary_path,
            output_tile_results,
            output_simple_summary_path,
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
    """Subtype microbial genomes using SNV targeting k-mer subtyping schemes.

    Includes subtyping schemes for:

    - Salmonella enterica spp. enterica serovar Heidelberg

    - Salmonella enterica spp. enterica serovar Enteritidis

    Developed by Geneviève Labbé, James Robertson, Peter Kruczkiewicz, Marisa Rankin, Matthew Gopez, Chad R. Laing, Philip Mabon, Kim Ziebell, Aleisha R. Reimer, Lorelee Tschetter, Gary Van Domselaar, Sadjia Bekal, Kimberley A. MacDonald, Linda Hoang, Linda Chui, Danielle Daignault, Durda Slavic, Frank Pollari, E. Jane Parmley, Philip Mabon, Elissa Giang, Lok Kan Lee, Jonathan Moffat, Marisa Rankin, Joanne MacKinnon, Roger Johnson, John H.E. Nash.
    """
    does_file_exist(output_simple_summary_path, force)
    does_file_exist(output_summary_path, force)
    does_file_exist(output_tile_results, force)
    scheme_fasta = get_scheme_fasta(scheme)
    scheme_subtype_counts = subtype_counts(scheme_fasta)

    subtyping_params = init_subtyping_params(**locals())
    input_contigs, input_reads = collect_inputs(**locals())
    if len(input_contigs) == 0 and len(input_reads) == 0:
        no_files_exception = click.UsageError('No input files specified!')
        click.secho('Please see -h/--help for more info', err=True)
        raise no_files_exception
    df_md = None
    if scheme_metadata:
        df_md = read_metadata_table(scheme_metadata)

    subtype_results = []  # type: List[Tuple[Subtype, pd.DataFrame]]
    if len(input_contigs) > 0:
        contigs_results = subtype_contigs_samples(input_genomes=input_contigs,
                                                  scheme=scheme,
                                                  scheme_name=scheme_name,
                                                  subtyping_params=subtyping_params,
                                                  scheme_subtype_counts=scheme_subtype_counts,
                                                  n_threads=threads)
        logging.info(f'Generated {len(contigs_results)} subtyping results from {len(input_contigs)} contigs samples')
        subtype_results += contigs_results
    if len(input_reads) > 0:
        reads_results = subtype_reads_samples(reads=input_reads,
                                              scheme=scheme,
                                              scheme_name=scheme_name,
                                              subtyping_params=subtyping_params,
                                              scheme_subtype_counts=scheme_subtype_counts,
                                              n_threads=threads)
        logging.info(f'Generated {len(reads_results)} subtyping results from {len(input_reads)} reads samples')
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
        logging.info(f'Wrote subtyping output summary to "{output_summary_path}"')
    else:
        # if no output path specified for the summary results, then print to stdout
        print(dfsummary.to_csv(sep='\t', index=None))

    if output_tile_results:
        if len(dfs) > 0:
            dfall = pd.concat(dfs)  # type: pd.DataFrame
            dfall.to_csv(output_tile_results, **kwargs_for_pd_to_table)
            logging.info(f'Tile results written to "{output_tile_results}".')
            if json_output:
                dfall.to_json(JSON_EXT_TMPL.format(output_tile_results), **kwargs_for_pd_to_json)
                logging.info(
                    f'Tile results written to "{JSON_EXT_TMPL.format(output_tile_results)}" in JSON format.')
        else:
            logging.error(
                f'No tile results generated. No tile results file written to "{output_tile_results}".')

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


def parse_comma_delimited_floats(ctx: click.Context, param: click.Option, value: Optional[str]) -> Optional[
    List[float]]:
    if value is None:
        return value

    if ',' in value:
        # TODO: more validation of user input values?
        return [float(x) for x in value.split(',')]
    else:
        return [float(value)]


@cli.command()
@click.option('-v', '--vcf-file-path',
              required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='Variant calling file (VCF) of a collection of input genomes for population of interest against a '
                   'reference genome that must be specified with --reference-genome-path')
@click.option('-r', '--reference-genome-path',
              required=True,
              type=click.Path(exists=True, dir_okay=False),
              help='Reference genome assembly file path. The reference used in the creation of the input VCF file.')
@click.option('--phylo-tree-path',
              required=False,
              type=click.Path(exists=True, dir_okay=False),
              help='Optional phylogenetic tree created from variant calling analysis.')
@click.option('-d', '--distance-thresholds',
              required=False,
              type=str,
              callback=parse_comma_delimited_floats,
              help='Comma delimited list of distance thresholds for creating hierarchical clustering groups '
                   '(e.g. "0,0.05,0.1,0.15")')
@click.option('-o', '--output-folder-name',
              required=True,
              type=click.Path(exists=False, dir_okay=True),
              help='Output folder name in which schema file would be located'
              )
@click.option('-s', '--schema-name',
              required=False,
              default="biohansel-schema",
              type=str,
              help='A unique name for the schema file that is generated, the default is just'
                   '{bio_hansel-schema}-reference_genome_name}-{schema_version}')
@click.option('-m', '--schema-version',
              required=False,
              default="0.1.0",
              type=str,
              help='An optional version number for the schema file that is generated')
@click.option('-t', '--tile-length',
              required=True,
              type=int,
              help='Length of sequence to be added around each SNV in the schema file,'
              'if an even integer is provided, then the tile length would be n+1'
              )
@click.option('-f', '--reference-genome-format',
              required=False,
              default=None,
              type=click.Choice(['fasta', 'genbank']),
              help='Reference genome file format: can be either fasta or genbank format'
              )
@click.option('-g', '--group-size-range',
              type=(int, int),
              default=(2,10),
              help='The range of child group size for each new subtype branching point from the parent group'
              )
@click.option('-p', '--pairwise-distance-metric',
                type=click.Choice(['hamming', 'euclidean', 'minkowski', 'cityblock', 'cosine', 'sqeuclidean', 'correlation', 'jaccard', 'chebyshev', 'braycurtis']),
                default='hamming',
                help='The distance metric used to calculate pairwise distances between SNVs'
                )
@click.option('-l', '--linkage-method',
                type=click.Choice(['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']),
                default='complete',
                help='The linkage method used to perform hierarchical clustering on SNVs'
                )

def create(vcf_file_path, reference_genome_path, phylo_tree_path, distance_thresholds, output_folder_name, schema_name,
           reference_genome_format, tile_length, group_size_range, schema_version, pairwise_distance_metric, linkage_method):
    """Create a biohansel subtyping scheme.

    From the results of a variant calling analysis, create a biohansel subtyping with single nucleotide variants (SNV)
    that discriminate subpopulations of genomes from all other genomes.
    """
    logging.info(f'VCF file path: {vcf_file_path}')
    logging.info(f'Reference genome file path: {reference_genome_path}')
    logging.info(f'Phylogenetic tree file path: {phylo_tree_path}')
    logging.info(f'Distance thresholds: {distance_thresholds}')
    logging.info(f'Output folder name: {output_folder_name}')
    logging.info(f'Scheme name: {schema_name}')
    logging.info(f'Reference genome format: {reference_genome_format}')
    logging.info(f'Group size range: {group_size_range}')
    logging.info(f'Padding Sequence length: {tile_length}')
    logging.info(f'Creating biohansel subtyping scheme from SNVs in "{vcf_file_path}" using reference genome ')
    logging.info(f'Pairwise distance metric: {pairwise_distance_metric}')
    logging.info(f'Linkage method to be used: {linkage_method}')
    reference_genome_name = genome_name_from_fasta_path(reference_genome_path)
    schema_name = f"{schema_name}-{reference_genome_name}-{schema_version}"

    if not os.path.exists(output_folder_name):
        os.makedirs(output_folder_name)
    

    sequence_df, binary_df = parse_vcf(vcf_file_path)
    logging.info(type(group_size_range))
    clusters = find_clusters(binary_df, group_size_range, distance_thresholds, pairwise_distance_metric, linkage_method)
   
    if phylo_tree_path is not None:
        phylo_tree_string=display_tree(phylo_tree_path, clusters.flat_clusters, output_folder_name)
    #sequence_records: Dict[str, Seq.Seq]
    sequence_records = parse_sequence_file(reference_genome_path, reference_genome_format)
    results_dict = group_snvs(binary_df, sequence_df, clusters.flat_clusters)
    for group, curr_df in results_dict.items():
        df_list = get_sequences(curr_df, tile_length,
                                sequence_records)
        write_sequence_file(output_folder_name, df_list, schema_name, group)
    output_schema_path = os.path.join(output_folder_name, f"{schema_name}.fasta")
    logging.info(f"Finished writing schema file to {output_schema_path}")
