# -*- coding: utf-8 -*-

from typing import Tuple, Optional

import pandas as pd

from ..qc.const import QC
from ..qc.utils import get_conflicting_tiles, get_num_pos_neg_tiles, get_mixed_subtype_tile_counts
from ..subtype import Subtype
from ..subtyping_params import SubtypingParams


def is_overall_coverage_low(st: Subtype, df: pd.DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    if not st.are_subtypes_consistent \
            or st.subtype is None \
            or not st.is_fastq_input():
        return None, None

    if st.avg_tile_coverage < p.min_coverage_warning:
        return QC.WARNING, f'Low coverage for all tiles ({st.avg_tile_coverage:.3f} < {p.min_coverage_warning} expected)'
    return None, None


def is_missing_tiles(st: Subtype, df: pd.DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """Are there more missing tiles than tolerated?

    Note:
        For reads, calculate the average coverage depth from the tiles that are present and provide an adequate
        error message based on the coverage.

    Args:
        st: Subtype results
        df: Subtyping results dataframe
        p: Subtyping/QC parameters specifically `max_perc_missing_tiles` for % missing tiles threshold

    Returns:
        None, None if less missing tiles than tolerate; otherwise, "FAIL", error message
    """

    if st.are_subtypes_consistent:
        return check_for_missing_tiles(is_fastq=st.is_fastq_input(),
                                       subtype_result=st.subtype,
                                       scheme=st.scheme,
                                       df=df,
                                       exp=int(st.n_tiles_matching_all_expected),
                                       obs=int(st.n_tiles_matching_all),
                                       p=p)
    else:
        message_list = []

        subtype_list = st.subtype.split(';')
        n_tiles_matching_expected = st.n_tiles_matching_all_expected.split(';')

        dfpos = df[df.is_pos_tile]
        mixed_subtype_counts = get_mixed_subtype_tile_counts(dfpos=dfpos,
                                                             subtype_list=subtype_list)

        tiles_matching_negative = st.n_tiles_matching_negative
        for curr_subtype, exp in zip(subtype_list, n_tiles_matching_expected):
            # We can omit the status because there will be a fail status already from non consistent subtypes.

            obs = mixed_subtype_counts.get(curr_subtype) + int(tiles_matching_negative)
            _, curr_messages = check_for_missing_tiles(is_fastq=st.is_fastq_input(),
                                                       subtype_result=curr_subtype,
                                                       scheme=st.scheme,
                                                       df=df,
                                                       exp=int(exp),
                                                       obs=obs,
                                                       p=p)

            message_list.append(curr_messages)

        error_messages = ' | '.join(filter(None.__ne__, message_list))

        return QC.FAIL, error_messages


def check_for_missing_tiles(is_fastq: bool,
                            subtype_result: str,
                            scheme: str,
                            df: pd.DataFrame,
                            exp: int,
                            obs: int,
                            p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """Check if there are too many missing tiles

    Also check if the mean tile coverage depth is above the low coverage threshold.

    Args:
        is_fastq: Is input sample reads?
        subtype_result: Single subtype designation
        scheme: Scheme name
        df: Subtyping results dataframe
        exp: Expected number of tiles that should be found
        obs: Actual observed number of tiles found
        p: Subtyping parameters

    Returns:
        Tuple of QC status and any QC messages
    """
    status = None
    messages = None

    # proportion of missing tiles
    p_missing = (exp - obs) / exp  # type: float
    if p_missing > p.max_perc_missing_tiles:
        status = QC.FAIL
        if is_fastq:
            tiles_with_hits = df[df['is_kmer_freq_okay']]  # type: pd.DataFrame
            depth = tiles_with_hits['freq'].mean()
            if depth < p.low_coverage_depth_freq:
                coverage_msg = f'Low coverage depth ({depth:.1f} < {float(p.low_coverage_depth_freq):.1f} expected); ' \
                               f'you may need more WGS data.'
            else:
                coverage_msg = f'Okay coverage depth ({depth:.1f} >= {float(p.low_coverage_depth_freq):.1f} expected), ' \
                               f'but this may be the wrong serovar or species for scheme "{scheme}"'
            messages = f'{p_missing:.2%} missing tiles; more than {p.max_perc_missing_tiles:.2%} missing ' \
                       f'tiles threshold. {coverage_msg}'
        else:
            messages = f'{p_missing:.2%} missing tiles for subtype "{subtype_result}"; more than ' \
                       f'{p.max_perc_missing_tiles:.2%} missing tile threshold'

    return status, messages


def is_mixed_subtype(st: Subtype, df: pd.DataFrame, *args) -> Tuple[Optional[str], Optional[str]]:
    """Is the subtype result mixed?

    Note:
        Check for `+` and `-` tiles that are present for the target site above minimum coverage threshold.
        If the subtyping result came back with an inconsistent call, provide an adequate error message to the user.

    Args:
        st: Subtype results
        df: Subtyping results dataframe
        args: unused args

    Returns:
        None, None if not mixed subtype result; otherwise, "FAIL", error message
    """
    if not st.are_subtypes_consistent:
        return QC.FAIL, f'Mixed subtypes found: "{"; ".join(sorted(st.inconsistent_subtypes))}".'
    conflicting_tiles = get_conflicting_tiles(st, df)
    if conflicting_tiles is None or conflicting_tiles.shape[0] == 0:
        return None, None

    s = 's' if conflicting_tiles.shape[0] > 1 else ''
    positions = ', '.join(conflicting_tiles['refposition'].astype(str).tolist())
    return QC.FAIL, f'Mixed subtype; the positive and negative tiles were found for ' \
                    f'the same target site{s} {positions} for subtype "{st.subtype}".'


def is_missing_too_many_target_sites(st: Subtype, df: pd.DataFrame, p: SubtypingParams) -> Tuple[
    Optional[str], Optional[str]]:
    """Are there too many missing target sites for an expected subtype?

    Check if there are any refpositions missing from the subtyping scheme in the result.
    If there are missing refpositions, this is a missing target.
    If there are too many missing targets, this result is an Ambiguous Result.

    Args:
        st: Subtype results
        df: Subtyping results dataframe
        p: Subtyping/QC parameters

    Returns:
        None, None if less missing targets than tolerated; otherwise, "FAIL", error message
    """
    if not st.are_subtypes_consistent \
            or st.subtype is None:
        return None, None

    potential_subtypes = st.all_subtypes.split('; ')
    uniq_positions = {y for x in potential_subtypes for y in st.scheme_subtype_counts[x].refpositions}
    missing_targets = [x for x in uniq_positions if (df.refposition == x).sum() == 0]

    exp = int(st.n_tiles_matching_all_expected)
    obs = int(st.n_tiles_matching_all)
    if (exp - obs) / exp <= p.max_perc_missing_tiles and len(missing_targets) >= p.min_ambiguous_tiles:
        return QC.FAIL, f'{QC.AMBIGUOUS_RESULTS_ERROR_3}: There were {len(missing_targets)} missing positions for ' \
                        f'subtype "{st.subtype}".'
    return None, None


def is_missing_downstream_targets(st: Subtype, *args) -> Tuple[Optional[str], Optional[str]]:
    """Are there missing downstream targets?

    Note:
        This method will check if there's any non_present_subtypes in the result, which would indicate a non confident
        result. This is due to the fact if you have a subtyping result of `2.1.1.2` and you're missing `2.1.1.2.X`
        You can't be sure that the subtype's final call is 2.1.1.2 if you're missing information.

    Args:
        st: Subtype results
        args: unused extra args

    Returns:
        None, None if no missing downstream targets; otherwise, "FAIL", error message
    """
    if st.non_present_subtypes:
        return QC.FAIL, f'{QC.UNCONFIDENT_RESULTS_ERROR_4}: Subtype "{st.subtype}" was found, but tiles for ' \
                        f'downstream subtype(s) "{st.non_present_subtypes}" were missing. Due to missing downstream ' \
                        f'tiles, there is a lack of confidence in the final subtype call.'
    return None, None


def is_maybe_intermediate_subtype(st: Subtype, df: pd.DataFrame, p: SubtypingParams) -> Tuple[
    Optional[str], Optional[str]]:
    """Is the result a possible intermediate subtype?

    Return a WARNING message if all the conditions are true:
    - 95% of the scheme tiles are found
    - 0 conflicting tiles
    - total subtype tiles < expected
    - +/- tiles exist in the result.

    Args:
        st: Subtype results
        df: Subtyping results dataframe
        p: Subtyping/QC parameters specifically using `max_perc_intermediate_tiles`

    Returns:
        None, None if no intermediate subtype possible; otherwise, "FAIL", error message
    """
    if not st.are_subtypes_consistent \
            or st.subtype is None:
        return None, None

    total_subtype_tiles = int(st.n_tiles_matching_subtype_expected)
    total_subtype_tiles_hits = int(st.n_tiles_matching_subtype)
    conflicting_tiles = get_conflicting_tiles(st, df)
    num_pos_tiles, num_neg_tiles = get_num_pos_neg_tiles(st, df)
    obs = int(st.n_tiles_matching_all)
    exp = int(st.n_tiles_matching_all_expected)
    if (exp - obs) / exp <= p.max_perc_intermediate_tiles and conflicting_tiles.shape[0] == 0 and \
            total_subtype_tiles_hits < total_subtype_tiles and num_pos_tiles and num_neg_tiles:
        return QC.WARNING, f'Possible intermediate subtype. All scheme tiles were found, but a fraction ' \
                           f'were positive for the final subtype. Total subtype matches observed ' \
                           f'(n={total_subtype_tiles_hits}) vs expected (n={total_subtype_tiles})'
    return None, None
