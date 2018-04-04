# -*- coding: utf-8 -*-

from typing import Tuple, Optional

from pandas import DataFrame

from ..subtyping_params import SubtypingParams
from ..qc.const import QC
from ..qc.utils import get_conflicting_tiles, get_num_pos_neg_tiles, get_mixed_subtype_tile_counts
from ..subtype import Subtype


def is_overall_coverage_low(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    if not st.are_subtypes_consistent \
            or st.subtype is None \
            or not st.is_fastq_input():
        return None, None

    if st.avg_tile_coverage < p.min_coverage_warning:
        return QC.WARNING, '{}: Low coverage for all tiles ({:.3f} < {} expected)'.format(
            QC.LOW_COVERAGE_WARNING,
            st.avg_tile_coverage,
            p.min_coverage_warning
        )

    return None, None


def is_missing_tiles(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
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
                                       curr_subtype=st.subtype,
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

        for curr_subtype, exp in zip(subtype_list, n_tiles_matching_expected):
            # We can omit the status because there will be a fail status already from non consistent subtypes.
            _, curr_messages = check_for_missing_tiles(is_fastq=st.is_fastq_input(),
                                                       curr_subtype=curr_subtype,
                                                       scheme=st.scheme,
                                                       df=df,
                                                       exp=int(exp),
                                                       obs=(mixed_subtype_counts.get(curr_subtype) +
                                                            int(st.n_tiles_matching_negative)),
                                                       p=p)

            message_list.append(curr_messages)

        error_messages = ' | '.join(filter(None.__ne__, message_list))

        return _, error_messages


def check_for_missing_tiles(is_fastq: bool, curr_subtype: str, scheme: str,
                            df: DataFrame, exp: int, obs: int, p: SubtypingParams):
    error_status = None
    error_messages = None

    p_missing = (exp - obs) / exp  # type: float
    if p_missing > p.max_perc_missing_tiles:
        if is_fastq:
            tiles_with_hits = df[df['is_kmer_freq_okay']]  # type: DataFrame
            depth = tiles_with_hits['freq'].mean()
            if depth < p.low_coverage_depth_freq:
                coverage_msg = 'Low coverage depth ({:.1f} < {:.1f} expected); you may need more WGS data.'.format(
                    depth,
                    float(p.low_coverage_depth_freq))
            else:
                coverage_msg = 'Okay coverage depth ({:.1f} >= {:.1f} expected), but this may be the wrong ' \
                               'serovar or species for scheme "{}"'.format(
                                depth,
                                float(p.low_coverage_depth_freq),
                                scheme)
            error_messages = '{status}: {p_missing:.2%} missing tiles; more than {p_missing_threshold:.2%} missing ' \
                             'tiles threshold. {coverage_msg}'.format(
                                status=QC.MISSING_TILES_ERROR_1,
                                p_missing=p_missing,
                                p_missing_threshold=p.max_perc_missing_tiles,
                                coverage_msg=coverage_msg)
        else:
            error_messages = '{}: {:.2%} missing tiles for subtype {}; more than {:.2%} missing tile threshold'.format(
                QC.MISSING_TILES_ERROR_1,
                p_missing,
                curr_subtype,
                p.max_perc_missing_tiles)
        error_status = QC.FAIL

    return error_status, error_messages


def is_mixed_subtype(st: Subtype, df: DataFrame, *args) -> Tuple[Optional[str], Optional[str]]:
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
        return QC.FAIL, \
               '{}: Mixed subtypes found: "{}".'.format(
                   QC.MIXED_SAMPLE_ERROR_2,
                   '; '.join(st.inconsistent_subtypes))
    conflicting_tiles = get_conflicting_tiles(st, df)
    if conflicting_tiles is None or conflicting_tiles.shape[0] == 0:
        return None, None

    return QC.FAIL, '{}: Mixed subtype detected. Positive and negative tiles detected for ' \
                    'the same target site "{}" for subtype "{}".'.format(
                        QC.MIXED_SAMPLE_ERROR_2,
                        '; '.join(conflicting_tiles['refposition'].astype(str).tolist()),
                        st.subtype)


def is_missing_too_many_target_sites(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[
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
        return QC.FAIL, \
               '{}: {} missing positions detected for subtype: {}.'.format(
                   QC.AMBIGUOUS_RESULTS_ERROR_3,
                   len(missing_targets),
                   st.subtype)
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
        return QC.FAIL, \
               '{}: Subtype {} was found, but tiles for downstream subtype(s) "{}" were missing. ' \
               'Due to missing downstream tiles, there is a lack of confidence in the final subtype call.'.format(
                   QC.UNCONFIDENT_RESULTS_ERROR_4,
                   st.subtype,
                   st.non_present_subtypes)

    return None, None


def is_maybe_intermediate_subtype(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[
    Optional[str], Optional[str]]:
    '''Is the result a possible intermediate subtype?

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
    '''
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
        return QC.WARNING, \
               '{}: Possible intermediate subtype. All scheme tiles were found, but a fraction ' \
               'were positive for the final subtype. Total subtype matches observed (n={}) vs expected (n={})'.format(
                   QC.INTERMEDIATE_SUBTYPE_WARNING,
                   total_subtype_tiles_hits,
                   total_subtype_tiles)
    return None, None
