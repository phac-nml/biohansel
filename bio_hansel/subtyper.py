# -*- coding: utf-8 -*-

import re
import logging
from typing import Optional, List, Dict, Union, Tuple

import attr
import pandas as pd

from .aho_corasick import init_automaton, find_in_fasta, find_in_fastqs
from .const import COLUMNS_TO_REMOVE
from .qc import perform_quality_check
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
    Subtype, Optional[pd.DataFrame]]:
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
        st.are_subtypes_consistent = False
        return st, None

    refpositions = [x for x, y in df.tilename.str.split('-')]
    subtypes = [y for x, y in df.tilename.str.split('-')]
    df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    df['subtype'] = subtypes
    df['is_pos_tile'] = ~df.tilename.str.contains('negative')
    dfpos = df[df.is_pos_tile]
    subtype_lens = dfpos.subtype.apply(len)
    max_subtype_strlen = subtype_lens.max()
    logging.debug('max substype str len: %s', max_subtype_strlen)
    dfpos_highest_res = dfpos[subtype_lens == max_subtype_strlen]
    logging.debug('dfpos_highest_res: %s', dfpos_highest_res)
    pos_subtypes = [[int(y) for y in x.split('.')] for x in dfpos.subtype.unique()]
    pos_subtypes.sort(key=lambda a: len(a))
    logging.debug('pos_subtypes: %s', pos_subtypes)
    inconsistent_subtypes = find_inconsistent_subtypes(pos_subtypes)
    logging.debug('inconsistent_subtypes: %s', inconsistent_subtypes)
    st.n_tiles_matching_all = df.tilename.unique().size
    st.n_tiles_matching_positive = dfpos.tilename.unique().size
    st.n_tiles_matching_subtype = dfpos_highest_res.tilename.unique().size
    pos_subtypes_str = [x for x in dfpos.subtype.unique()]
    pos_subtypes_str.sort(key=lambda x: len(x))
    st.all_subtypes = '; '.join(pos_subtypes_str)
    subtype_list = [x for x in dfpos_highest_res.subtype.unique()]
    st.subtype = '; '.join(subtype_list)
    st.n_tiles_matching_all_expected = ';'.join([str(scheme_subtype_counts[x].all_tile_count) for x in subtype_list])
    st.n_tiles_matching_positive_expected = ';'.join(
        [str(scheme_subtype_counts[x].positive_tile_count) for x in subtype_list])
    st.n_tiles_matching_subtype_expected = ';'.join(
        [str(scheme_subtype_counts[x].subtype_tile_count) for x in subtype_list])
    st.tiles_matching_subtype = '; '.join([x for x in dfpos_highest_res.tilename.unique()])

    if len(inconsistent_subtypes) > 0:
        st.are_subtypes_consistent = False
        st.inconsistent_subtypes = inconsistent_subtypes

    possible_downstream_subtypes = [s for s in scheme_subtype_counts
                                    if re.search("^({})(\.)(\d)$".format(re.escape(st.subtype)), s)]
    st.non_present_subtypes = [x for x in possible_downstream_subtypes
                               if not (df.subtype == x).any()]
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
        subtype_results.append(attr.asdict(subtype))


def parallel_query_reads_ac(reads: List[Tuple[List[str], str]],
                            scheme: str,
                            scheme_name: Optional[str] = None,
                            subtyping_params: Optional[SubtypingParams] = None,
                            scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                            n_threads: int = 1) -> List[Tuple[Subtype, pd.DataFrame]]:
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
    logging.info("genome_name %s", genome_name)

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
        return st, None

    refpositions = [x for x, y in df.tilename.str.split('-')]
    subtypes = [y for x, y in df.tilename.str.split('-')]
    df['refposition'] = [int(x.replace('negative', '')) for x in refpositions]
    df['subtype'] = subtypes
    df['is_pos_tile'] = ~df.tilename.str.contains('negative')
    df['is_kmer_freq_okay'] = (df.freq >= subtyping_params.min_kmer_freq) & (df.freq <= subtyping_params.max_kmer_freq)
    dfgood = df[df.is_kmer_freq_okay]
    dfpos = dfgood[dfgood.is_pos_tile]
    dfneg = dfgood[~dfgood.is_pos_tile]
    subtype_lens = dfpos.subtype.apply(len)
    max_subtype_strlen = subtype_lens.max()
    logging.debug('max substype str len: %s', max_subtype_strlen)
    dfpos_highest_res = dfpos[subtype_lens == max_subtype_strlen]
    pos_subtypes = [[int(y) for y in x.split('.')] for x in dfpos.subtype.unique()]
    pos_subtypes.sort(key=lambda a: len(a))
    logging.debug('pos_subtypes: %s', pos_subtypes)
    inconsistent_subtypes = find_inconsistent_subtypes(pos_subtypes)
    logging.debug('inconsistent_subtypes: %s', inconsistent_subtypes)
    st.n_tiles_matching_all = dfgood.shape[0]
    st.n_tiles_matching_positive = dfpos.shape[0]
    st.n_tiles_matching_negative = dfneg.shape[0]
    st.n_tiles_matching_subtype = dfpos_highest_res.shape[0]
    pos_subtypes_str = [x for x in dfpos.subtype.unique()]
    pos_subtypes_str.sort(key=lambda x: len(x))
    st.all_subtypes = '; '.join(pos_subtypes_str)
    subtype_list = [x for x in dfpos_highest_res.subtype.unique()]
    st.subtype = '; '.join(subtype_list)
    st.n_tiles_matching_all_expected = ';'.join([str(scheme_subtype_counts[x].all_tile_count) for x in subtype_list])
    st.n_tiles_matching_positive_expected = ';'.join(
        [str(scheme_subtype_counts[x].positive_tile_count) for x in subtype_list])
    st.n_tiles_matching_subtype_expected = ';'.join(
        [str(scheme_subtype_counts[x].subtype_tile_count) for x in subtype_list])
    st.tiles_matching_subtype = '; '.join([x for x in dfpos_highest_res.tilename.unique()])
    if len(inconsistent_subtypes) > 0:
        st.are_subtypes_consistent = False
        st.inconsistent_subtypes = inconsistent_subtypes
    possible_downstream_subtypes = [s for s in scheme_subtype_counts
                                    if re.search("^({})(\.)(\d)$".format(re.escape(st.subtype)), s)]
    st.non_present_subtypes = [x for x in possible_downstream_subtypes
                               if not (df.subtype == x).any()]
    st.avg_tile_coverage = df['freq'].mean()
    st.qc_status, st.qc_message = perform_quality_check(st, df, subtyping_params)
    df['sample'] = genome_name
    df['scheme'] = scheme_name or scheme
    df['scheme_version'] = scheme_version
    df['qc_status'] = st.qc_status
    df['qc_message'] = st.qc_message
    df = df[df.columns[~df.columns.isin(COLUMNS_TO_REMOVE)]]
    return st, df


def query_reads_ac(subtype_results: List[Subtype],
                   dfs: List[pd.DataFrame],
                   reads: List[Tuple[List[str], str]],
                   scheme: str,
                   scheme_name: Optional[str] = None,
                   subtyping_params: Optional[SubtypingParams] = None,
                   scheme_subtype_counts: Optional[Dict[str, SubtypeCounts]] = None,
                   n_threads: int = 1):
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
