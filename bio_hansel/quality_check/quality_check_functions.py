from pandas import DataFrame

from ..subtyping_params import SubtypingParams
from ..quality_check.const import FAIL_MESSAGE, WARNING_MESSAGE, MISSING_TILES_ERROR_1, \
    MIXED_SAMPLE_ERROR_2, AMBIGUOUS_RESULTS_ERROR_3, INTERMEDIATE_SUBTYPE_WARNING, \
    NON_CONFIDENT_RESULTS_ERROR_4
from ..quality_check.qc_utils import get_conflicting_tiles, get_num_pos_neg_tiles
from ..subtype import Subtype
from typing import Tuple, Optional
import logging


def does_subtype_result_exist(st: Subtype) -> bool:
    """ Check if the subtype result exists.
    Note:
            This method verifies that there is an expected matching subtype so we can proceed
            with quality checking.

    Input:
            :param st: Subtyping results.

    Returns:
            Bool: True if the subtype exists
            Bool: False if the subtype does not exist
    """
    logging.debug("QC: Checking if subtype result exists.")
    return st.subtype is not None and 0 < len(st.subtype)


def check_missing_tiles(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """ Check if there's more than 5% missing tiles.
    Note:
            check_missing_tiles will check if more than 5% of the scheme's tiles are missing.
            If they are, we calculate the average frequency coverage depth for the tiles that are present, and
            provide an adequate error message based on the coverage.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.
            :param p: Subtyping parameters for quality check thresholds.

    Returns:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning, Fail }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for missing tiles.")
    error_status = None
    error_messages = None

    if not st.are_subtypes_consistent:
        logging.debug("QC: Checking for missing tiles not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for missing tiles not run."
        error_status = FAIL_MESSAGE
        return error_status, error_messages

    exp = int(st.n_tiles_matching_all_expected)
    obs = int(st.n_tiles_matching_all)

    if (exp - obs) / exp > p.max_perc_missing_tiles:
        tiles_with_hits = df[df['is_kmer_freq_okay']] # type: DataFrame
        avg_depth = tiles_with_hits['freq'].mean()

        error_messages = "{}: More than {:.2%} missing tiles were detected. {} Avg calculated tile coverage = {}".format(
            MISSING_TILES_ERROR_1, p.max_perc_missing_tiles,
            "Low coverage detected, possibly need more whole genome sequencing data."
            if avg_depth < p.low_coverage_depth_freq
            else "Adequate coverage detected, this may be the wrong serovar/species for scheme: {}".format(st.scheme),
            avg_depth
        )

        error_status = FAIL_MESSAGE

    return error_status, error_messages


def check_mixed_subtype(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """ Check to see if the subtype is mixed.
    Note:
            check_mixed_subtype will check the results of the typing for:
            1) + and - tiles that are present for the target site above minimum coverage threshold.
            2) If the subtyping result came back with an inconsistent call.
            Then the method will provide an adequate error message to the user.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.
            :param p: Subtyping parameters for quality check thresholds.

    Returns:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning, Fail }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for mixed subtypes.")
    error_status = None
    error_messages = None

    # First check the obvious, if the st is consistent, of course it's a mixed subtype.
    if not st.are_subtypes_consistent or st.inconsistent_subtypes:
        error_status = FAIL_MESSAGE
        error_messages = "{}: Mixed subtypes detected. "\
                         "Mixed subtypes found: {}.".format(MIXED_SAMPLE_ERROR_2, st.inconsistent_subtypes)
    else:
        # Find the conflicting tiles
        conflicting_tiles = get_conflicting_tiles(st, df)
        if 0 < conflicting_tiles.shape[0]:
            error_status = FAIL_MESSAGE
            error_messages = "{}: Mixed subtype detected." \
                             " Positive and negative tiles detected for the same target site" \
                             " {} for subtype {}.".format(MIXED_SAMPLE_ERROR_2,
                                                          conflicting_tiles['refposition'].tolist(), st.subtype)

    return error_status, error_messages


def is_missing_target_sites(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """Are there X or more missing target sites for an expected subtype? (ERROR 3)
    Note:
            This method will check if there are any refpositions missing from the subtyping scheme in the result.
            If there are missing refpositions, this is a missing target. If there are too many missing targets,
            this result is an Ambiguous Result.
    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.
            :param p: Subtyping parameters for quality check thresholds.

    Returns:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning, Fail }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for ambiguous results.")
    error_status = None
    error_messages = None

    # We do this check as inconsistent subtypes fall within mixed subtypes. This checks for inconsistent results.
    if not st.are_subtypes_consistent:
        logging.debug("QC: Checking for ambiguous results not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for ambiguous results not run."
        error_status = FAIL_MESSAGE
        return error_status, error_messages

    potential_subtypes = st.all_subtypes.split('; ')
    uniq_positions = {y for x in potential_subtypes for y in st.scheme_subtype_counts[x].refpositions}
    missing_targets = [x for x in uniq_positions if (df.refposition == x).sum() == 0]

    exp = int(st.n_tiles_matching_all_expected)
    obs = int(st.n_tiles_matching_all)
    if (exp - obs) / exp <= p.max_perc_missing_tiles and len(missing_targets) >= p.min_ambiguous_tiles:
            error_status = FAIL_MESSAGE
            error_messages = ("{}: {} missing positions detected for subtype: {}."
                              .format(
                                        AMBIGUOUS_RESULTS_ERROR_3,
                                        len(missing_targets),
                                        st.subtype
                                     )
                              )

    return error_status, error_messages


def is_missing_downstream_targets(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """ Check for missing downstream targets (ERROR 4)
    Note:
            This method will check if there's any non_present_subtypes in the result, which would indicate a non confident
            result. This is due to the fact if you have a subtyping result of `2.1.1.2` and you're missing `2.1.1.2.X`
            You can't be sure that the subtype's final call is 2.1.1.2 if you're missing information.

    Args:
            :param st: Subtyping results.

    Output:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { FAIL_MESSAGE }
            error_messages: Contains the reasoning for any fail status.
    """
    error_status = None
    error_messages = None

    non_present_subtypes = st.non_present_subtypes
    if non_present_subtypes:
        error_status = FAIL_MESSAGE
        error_messages = ("{}: Subtype {} was found, but downstream subtype(s) {} tiles were missing."
                          " Therefore due to the missing downstream tiles, there is a lack of confidence in the"
                          "final subtype call."
                          .format(
                            NON_CONFIDENT_RESULTS_ERROR_4,
                            st.subtype,
                            non_present_subtypes
                           )
                          )

    return error_status, error_messages


def check_intermediate_subtype(st: Subtype, df: DataFrame, p: SubtypingParams) -> Tuple[Optional[str], Optional[str]]:
    """ Check for intermediate subtypes within the result.
    Note:
            Check_intermediate_subtype will return a WARNING message if 95% of the scheme tiles are found, there are
            0 contradicting tiles (tiles with same position but there exists + and -), total subtype tiles < expected,
            and +/- tiles exist in the result.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.
            :param p: Subtyping parameters for quality check thresholds.

    Output:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning }
            error_messages: Contains the reasoning for any warning status.
    """
    logging.debug("QC: Checking for intermediate subtypes.")
    error_status = None
    error_messages = None

    if not st.are_subtypes_consistent:
        logging.debug("QC: Checking for intermediate subtypes not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for intermediate subtype not run."
        error_status = FAIL_MESSAGE
        return error_status, error_messages

    total_subtype_tiles = int(st.n_tiles_matching_subtype_expected)
    total_subtype_tiles_hits = int(st.n_tiles_matching_subtype)
    conflicting_tiles = get_conflicting_tiles(st, df)
    num_pos_tiles, num_neg_tiles = get_num_pos_neg_tiles(st, df)

    obs = int(st.n_tiles_matching_all)
    exp = int(st.n_tiles_matching_all_expected)

    if (exp - obs) / exp <= p.max_perc_intermediate_tiles and conflicting_tiles.shape[0] == 0 and \
            total_subtype_tiles_hits < total_subtype_tiles and num_pos_tiles and num_neg_tiles:

        error_status = WARNING_MESSAGE
        error_messages = "{}: Possible intermediate subtype. All scheme tiles were found, " \
                         "but a fraction were positive for the final subtype. " \
                         "Total subtype hits: {} | Total subtype expected: {}." \
                         .format(INTERMEDIATE_SUBTYPE_WARNING, total_subtype_tiles_hits, total_subtype_tiles)

    return error_status, error_messages
