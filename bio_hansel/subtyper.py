# -*- coding: utf-8 -*-
"""
Functions for subtyping of reads (e.g. FASTQ) and contigs (e.g. FASTA) using bio_hansel-compatible subtyping schemes.
"""
import re
import logging
from typing import Optional, List, Dict, Union, Tuple

import attr
import pandas as pd

from .aho_corasick import init_automaton, find_in_fasta, find_in_fastqs
from .const import COLUMNS_TO_REMOVE
from .qc import perform_quality_check, QC
from .subtype import Subtype
from .subtype_stats import SubtypeCounts
from .subtype_stats import subtype_counts
from .subtyping_params import SubtypingParams
from .utils import find_inconsistent_subtypes, get_scheme_fasta, get_scheme_version, init_subtyping_params


def subtype_contigs_ac(fasta_path: str,
                       genome_name: str,
                       scheme: str,
                       subtyping_params: Optional[SubtypingParams] = None,
                       scheme_name: Optional[str] = None,
                       scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None) -> Tuple[
    Subtype, pd.DataFrame]:
    """Subtype input contigs using a particular scheme.

    Args:
        fasta_path: Input FASTA file path
        genome_name: Input genome name
        scheme: bio_hansel scheme FASTA path
        subtyping_params: scheme specific subtyping parameters
        scheme_name: optional scheme name
        scheme_subtype_counts: summary information about scheme

    Returns:
        - Subtype result
        - pd.DataFrame of detailed subtyping results
    """
    scheme_fasta = get_scheme_fasta(scheme)
    if scheme_subtype_counts is None:
        scheme_subtype_counts = subtype_counts(scheme_fasta)
    if subtyping_params is None:
        subtyping_params = init_subtyping_params(scheme=scheme)
    scheme_version = get_scheme_version(scheme)
    st = Subtype(sample=genome_name,
                 file_path=fasta_path,
                 scheme=scheme_name or scheme,
                 scheme_version=scheme_version,
                 scheme_subtype_counts=scheme_subtype_counts)

    A = init_automaton(scheme_fasta)
    df = find_in_fasta(A, fasta_path)

    if df is None or df.shape[0] == 0:
        logging.warning('No subtyping tile matches for input "%s" for scheme "%s"', fasta_path, scheme)
        st.qc_status = QC.FAIL
        st.qc_message = QC.NO_SUBTYPE_RESULT
        st.are_subtypes_consistent = False
        return st, empty_results(st)

    refpositions = [x for x, y in df.tilename.str.split('-')]
    subtypes = [y for x, y in df.tilename.str.split('-')]
    df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    df['subtype'] = subtypes
    df['is_pos_tile'] = ~df.tilename.str.contains('negative')
    process_contigs_subtyping_results(st, df, scheme_subtype_counts)
    st.qc_status, st.qc_message = perform_quality_check(st, df, subtyping_params)

    logging.info(st)

    df['sample'] = genome_name
    df['file_path'] = fasta_path
    df['scheme'] = scheme_name or scheme
    df['scheme_version'] = scheme_version
    df['qc_status'] = st.qc_status
    df['qc_message'] = st.qc_message

    df = df[df.columns[~df.columns.isin(COLUMNS_TO_REMOVE)]]
    return st, df




def parallel_query_contigs_ac(input_genomes: List[Tuple[str, str]],
                              scheme: str,
                              scheme_name: Optional[str] = None,
                              subtyping_params: Optional[SubtypingParams] = None,
                              scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                              n_threads: int = 1):
    """Parallel subtyping of input contigs

    Subtype and analyse each input in parallel using a multiprocessing thread pool.

    Args:
        input_genomes: Input genome FASTA paths
        scheme: bio_hansel scheme FASTA path to subtype with
        scheme_name: optional scheme name
        subtyping_params: scheme specific subtyping parameters
        scheme_subtype_counts: scheme summary information
        n_threads: Number of threads to use

    Returns:
        A list of tuples of Subtype results and a pd.DataFrame of detailed subtyping results for each input
    """
    from multiprocessing import Pool
    logging.info('Initializing thread pool with %s threads', n_threads)
    pool = Pool(processes=n_threads)
    logging.info('Running analysis asynchronously on %s input genomes', len(input_genomes))
    res = [pool.apply_async(subtype_contigs_ac, (input_fasta,
                                                 genome_name,
                                                 scheme,
                                                 subtyping_params,
                                                 scheme_name,
                                                 scheme_subtype_counts))
           for input_fasta, genome_name in input_genomes]
    logging.info('Parallel analysis complete! Retrieving analysis results')
    outputs = [x.get() for x in res]
    return outputs


def query_contigs_ac(subtype_results: List[Subtype],
                     dfs: List[pd.DataFrame],
                     input_genomes: List[Tuple[str, str]],
                     scheme: str,
                     scheme_name: Optional[str] = None,
                     subtyping_params: Optional[SubtypingParams] = None,
                     scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                     n_threads: int = 1) -> None:
    """Subtype input genomes using a scheme.

    Args:
        subtype_results: Subtyping results for each input genome
        dfs: detailed subtyping results for each input genome
        input_genomes: input genomes; tuple of FASTA file path and genome name
        scheme: bio_hansel scheme FASTA path
        scheme_name: optional scheme name
        subtyping_params: scheme specific subtyping parameters
        scheme_subtype_counts: summary information about scheme
        n_threads: number of threads to use for subtyping analysis

    Returns:
        `subtype_results` and `dfs` by reference
    """
    if n_threads == 1:
        logging.info('Serial single threaded run mode on %s input genomes', len(input_genomes))
        outputs = [subtype_contigs_ac(fasta_path=input_fasta,
                                      genome_name=genome_name,
                                      scheme=scheme,
                                      scheme_name=scheme_name,
                                      subtyping_params=subtyping_params,
                                      scheme_subtype_counts=scheme_subtype_counts)
                   for input_fasta, genome_name in input_genomes]
    else:
        outputs = parallel_query_contigs_ac(input_genomes, scheme, scheme_name, subtyping_params, scheme_subtype_counts,
                                            n_threads)
    for subtype, df in outputs:
        if df is not None:
            dfs.append(df)
        else:
            logging.error(f'Subtyping results DataFrame is "None" for {subtype}. Building empty results DF.')
            dfs.append(empty_results(subtype))
        subtype_results.append(attr.asdict(subtype))


def parallel_query_reads_ac(reads: List[Tuple[List[str], str]],
                            scheme: str,
                            scheme_name: Optional[str] = None,
                            subtyping_params: Optional[SubtypingParams] = None,
                            scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                            n_threads: int = 1) -> List[Tuple[Subtype, pd.DataFrame]]:
    """Parallel subtyping of input reads

    Subtype and analyse each input in parallel using a multiprocessing thread pool.

    Args:
        reads: Input reads; list of tuples of FASTQ file paths and genome names
        scheme: bio_hansel scheme FASTA path
        scheme_name: optional scheme name
        subtyping_params: scheme specific subtyping parameters
        scheme_subtype_counts: scheme summary information
        n_threads: number of threads to use

    Returns:
        A list of tuples of Subtype results and a pd.DataFrame of detailed subtyping results for each input
    """
    from multiprocessing import Pool
    logging.info('Initializing thread pool with %s threads', n_threads)
    pool = Pool(processes=n_threads)
    logging.info('Running analysis asynchronously on %s input genomes', len(reads))
    res = [pool.apply_async(subtype_reads_ac, (fastqs,
                                               genome_name,
                                               scheme,
                                               scheme_name,
                                               subtyping_params,
                                               scheme_subtype_counts))
           for fastqs, genome_name in reads]
    logging.info('Parallel analysis complete! Retrieving analysis results')
    outputs = [x.get() for x in res]
    return outputs


def subtype_reads_ac(reads: Union[str, List[str]],
                     genome_name: str,
                     scheme: str,
                     scheme_name: Optional[str] = None,
                     subtyping_params: Optional[SubtypingParams] = None,
                     scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None) -> Tuple[
    Subtype, Optional[pd.DataFrame]]:
    """Subtype input reads using a particular scheme.

    Args:
        reads: Input FASTQ file path(s)
        genome_name: Input genome name
        scheme: bio_hansel scheme FASTA path
        scheme_name: optional scheme name
        subtyping_params: scheme specific subtyping parameters
        scheme_subtype_counts: summary information about scheme

    Returns:
        - Subtype result
        - pd.DataFrame of detailed subtyping results
    """
    scheme_fasta = get_scheme_fasta(scheme)
    if scheme_subtype_counts is None:
        scheme_subtype_counts = subtype_counts(scheme_fasta)
    if subtyping_params is None:
        subtyping_params = init_subtyping_params(scheme=scheme)
    scheme_version = get_scheme_version(scheme)

    st = Subtype(sample=genome_name,
                 file_path=reads,
                 scheme=scheme_name or scheme,
                 scheme_version=scheme_version,
                 scheme_subtype_counts=scheme_subtype_counts)

    A = init_automaton(scheme_fasta)
    if isinstance(reads, str):
        df = find_in_fastqs(A, reads)
    elif isinstance(reads, list):
        df = find_in_fastqs(A, *reads)
    else:
        raise ValueError('Unexpected type "{}" for "reads": {}'.format(type(reads), reads))

    if df is None or df.shape[0] == 0:
        logging.warning('No subtyping tile matches for input "%s" for scheme "%s"', reads, scheme)
        st.are_subtypes_consistent = False
        st.qc_status = QC.FAIL
        st.qc_message = QC.NO_SUBTYPE_RESULT
        return st, empty_results(st)

    refpositions = [x for x, y in df.tilename.str.split('-')]
    subtypes = [y for x, y in df.tilename.str.split('-')]
    df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    df['subtype'] = subtypes
    df['is_pos_tile'] = ~df.tilename.str.contains('negative')
    df['is_kmer_freq_okay'] = (df.freq >= subtyping_params.min_kmer_freq) & (df.freq <= subtyping_params.max_kmer_freq)
    st.avg_tile_coverage = df['freq'].mean()
    process_contigs_subtyping_results(st, df[df.is_kmer_freq_okay], scheme_subtype_counts)
    st.qc_status, st.qc_message = perform_quality_check(st, df, subtyping_params)
    df['sample'] = genome_name
    df['scheme'] = scheme_name or scheme
    df['scheme_version'] = scheme_version
    df['qc_status'] = st.qc_status
    df['qc_message'] = st.qc_message
    df = df[df.columns[~df.columns.isin(COLUMNS_TO_REMOVE)]]
    return st, df


def process_contigs_subtyping_results(st, df, scheme_subtype_counts):
    dfpos = df[df.is_pos_tile]
    dfpos_highest_res = highest_resolution_subtype_results(dfpos)
    subtype_list = [x for x in dfpos_highest_res.subtype.unique()]
    set_subtype_results(st, dfpos, subtype_list)
    set_inconsistent_subtypes(st, find_inconsistent_subtypes(sorted_positive_subtypes(dfpos)))
    calculate_subtyping_stats(st, df, dfpos, subtype_list, scheme_subtype_counts)
    set_non_present_subtypes(st, df, scheme_subtype_counts)


def set_subtype_results(st, dfpos, subtype_list):
    st.subtype = '; '.join(subtype_list)
    st.tiles_matching_subtype = '; '.join(subtype_list)
    pos_subtypes_str = [x for x in dfpos.subtype.unique()]
    pos_subtypes_str.sort(key=lambda x: len(x))
    st.all_subtypes = '; '.join(pos_subtypes_str)


def calculate_subtyping_stats(st, df, dfpos, subtype_list, scheme_subtype_counts):
    st.n_tiles_matching_all = df.tilename.unique().size
    st.n_tiles_matching_positive = dfpos.tilename.unique().size
    st.n_tiles_matching_negative = df[~df.is_pos_tile]
    st.n_tiles_matching_subtype = len(subtype_list)
    st.n_tiles_matching_all_expected = ';'.join([str(scheme_subtype_counts[x].all_tile_count) for x in subtype_list])
    st.n_tiles_matching_positive_expected = ';'.join(
        [str(scheme_subtype_counts[x].positive_tile_count) for x in subtype_list])
    st.n_tiles_matching_subtype_expected = ';'.join(
        [str(scheme_subtype_counts[x].subtype_tile_count) for x in subtype_list])


def count_periods(subtype_with_periods: str) -> pd.Series:
    return (pd.Series(list(subtype_with_periods)) == '.').sum()


def highest_resolution_subtype_results(dfpos):
    subtype_lens = dfpos.subtype.apply(count_periods)
    max_subtype_strlen = subtype_lens.max()
    logging.debug('max substype str len: %s', max_subtype_strlen)
    dfpos_highest_res = dfpos[subtype_lens == max_subtype_strlen]
    return dfpos_highest_res


def sorted_positive_subtypes(dfpos: pd.DataFrame) -> List[List[int]]:
    pos_subtypes = [[int(y) for y in x.split('.')] for x in dfpos.subtype.unique()]
    pos_subtypes.sort(key=lambda a: len(a))
    return pos_subtypes


def set_non_present_subtypes(st, df, scheme_subtype_counts):
    possible_downstream_subtypes = [s for s in scheme_subtype_counts
                                    if re.search("^({})(\.)(\d)$".format(re.escape(st.subtype)), s)]
    st.non_present_subtypes = [x for x in possible_downstream_subtypes
                               if not (df.subtype == x).any()]


def set_inconsistent_subtypes(st: Subtype, inconsistent_subtypes: List):
    if len(inconsistent_subtypes) > 0:
        st.are_subtypes_consistent = False
        st.inconsistent_subtypes = inconsistent_subtypes
    else:
        st.are_subtypes_consistent = True
        st.are_subtypes_consistent = None


def query_reads_ac(subtype_results: List[Subtype],
                   dfs: List[pd.DataFrame],
                   reads: List[Tuple[List[str], str]],
                   scheme: str,
                   scheme_name: Optional[str] = None,
                   subtyping_params: Optional[SubtypingParams] = None,
                   scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                   n_threads: int = 1):
    """Subtype input genomes using a scheme.

    Args:
        subtype_results: Subtyping results for each input genome
        dfs: detailed subtyping results for each input genome
        input_genomes: input genomes; tuple of FASTA file path and genome name
        scheme: bio_hansel scheme FASTA path
        scheme_name: optional scheme name
        subtyping_params: scheme specific subtyping parameters
        scheme_subtype_counts: summary information about scheme
        n_threads: number of threads to use for subtyping analysis

    Returns:
        `subtype_results` and `dfs` by reference
    """
    if n_threads == 1:
        logging.info('Serial single threaded run mode on %s input genomes', len(reads))
        outputs = [subtype_reads_ac(reads=fastq_files,
                                    genome_name=genome_name,
                                    scheme=scheme,
                                    scheme_name=scheme_name,
                                    subtyping_params=subtyping_params,
                                    scheme_subtype_counts=scheme_subtype_counts)
                   for fastq_files, genome_name in reads]
    else:
        outputs = parallel_query_reads_ac(reads=reads,
                                          scheme=scheme,
                                          scheme_name=scheme_name,
                                          subtyping_params=subtyping_params,
                                          scheme_subtype_counts=scheme_subtype_counts,
                                          n_threads=n_threads)
    for subtype, df in outputs:
        if df is not None:
            dfs.append(df)
        subtype_results.append(attr.asdict(subtype))


def empty_results(st: Subtype) -> pd.DataFrame:
    """Get a results DataFrame with 1 entry when no results are retrieved.

    When the output DataFrame for subtyping is empty, generate a "dummy" DataFrame with one row showing that subtyping
    was performed on the sample, but no results were generated because nothing was found.

    Args:
        st: Subtype results; should be a null result

    Returns:
        pd.DataFrame with 1 row with null result
    """
    return pd.DataFrame({0: dict(sample=st.sample,
                                 file_path=st.file_path,
                                 subtype=st.subtype,
                                 refposition=None,
                                 is_pos_tile=None,
                                 scheme=st.scheme,
                                 scheme_version=st.scheme_version,
                                 qc_status=st.qc_status,
                                 qc_message=st.qc_message)}).transpose()