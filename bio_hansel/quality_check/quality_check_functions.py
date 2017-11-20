from pandas import DataFrame
from bio_hansel.quality_check.const import FAIL_MESSAGE, WARNING_MESSAGE, MIN_TILES_THRESHOLD, MISSING_TILES_ERROR_1A, \
    MISSING_TILES_ERROR_1B, MIXED_SAMPLE_ERROR_2A, INCONSISTENT_RESULTS_ERROR_3A, INTERMEDIATE_SUBTYPE_WARNING, \
    INCONSISTENT_RESULTS_ERROR_3B
from bio_hansel.quality_check.qc_utils import get_conflicting_tiles, get_num_pos_neg_tiles, possible_subtypes_exist_in_df
from bio_hansel.subtype import Subtype
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
    return st.subtype is not None and len(st.subtype) > 0


def check_missing_tiles(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """ Check if there's more than 5% missing tiles.
    Note:   check_missing_tiles will check if more than 5% of the scheme's tiles are missing.
            If they are, we calculate the average frequency coverage depth for the tiles that are present, and
            provide an adequate error message based on the coverage.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning, Fail }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for missing tiles.")
    error_status = None
    error_messages = None

    if st.are_subtypes_consistent:
        total_tiles = int(st.n_tiles_matching_all_expected)
        total_tiles_hits = int(st.n_tiles_matching_all)

        if total_tiles_hits < (total_tiles - total_tiles * MIN_TILES_THRESHOLD):
            tiles_with_hits = df[(df['is_kmer_freq_okay'] == True)]
            average_freq_coverage_depth = tiles_with_hits['freq'].sum() / len(tiles_with_hits)

            # Per scheme error tailored messages.
            if str(st.scheme).lower() == 'heidelberg':
                if average_freq_coverage_depth < 20:
                    error_messages = " {}: More than 5% missing tiles were detected." \
                                     " Low coverage detected, possibly need more whole genome sequencing data." \
                                     " Average calculated tile coverage = {}"\
                                     .format(MISSING_TILES_ERROR_1A, str(average_freq_coverage_depth))
                else:
                    error_messages = "{}: More than 5% missing tiles were detected. " \
                                     "Adequate coverage detected, this may be the wrong serovar/species for scheme: {}"\
                                     " Average calculated tile coverage = {}"\
                                     .format(MISSING_TILES_ERROR_1B, str(st.scheme), average_freq_coverage_depth)
            else:
                error_messages = "More than 5% missing tiles were detected." \
                                 " Average calculated tile coverage = {}".format(str(average_freq_coverage_depth))

            error_status = FAIL_MESSAGE
    else:
        logging.debug("QC: Checking for missing tiles not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for missing tiles not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages


def check_mixed_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """ Check to see if the subtype is mixed.
    Note:
            check_mixed_subtype will check the results of the typing for:
            1) + and - tiles that are present for the target site above minimum coverage threshold.
            2) If the subtyping result came back with an inconsistent call.
            Then the method will provide an adequate error message to the user.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

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
                         "Mixed subtypes found: {}.".format(MIXED_SAMPLE_ERROR_2A, st.inconsistent_subtypes)
    else:
        # Find the conflicting tiles
        conflicting_tiles = get_conflicting_tiles(st, df)
        if conflicting_tiles:
            error_status = FAIL_MESSAGE
            error_messages = "{}: Mixed subtype detected." \
                             " Positive and negative tiles detected for the same target site" \
                             " {} for subtype {}.".format(MIXED_SAMPLE_ERROR_2A, conflicting_tiles, st.subtype)

    return error_status, error_messages


def check_inconsistent_results(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """ Check if the subtyping results are inconsistent results.
    Note:
            An inconsistent results error occurs when 3 or more tiles are missing (+ / -) for the subtype call,
            and there must be <= 5% tiles missing from the subtyping scheme.
            This method returns the error_status == FAIL if inconsistent results were found,
            and the error_messages corresponding to the error.
    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning, Fail }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for inconsistent results.")
    error_status = None
    error_messages = None

    # We do this check as inconsistent subtypes fall within mixed subtypes. This checks for inconsistent results.
    if st.are_subtypes_consistent:

        tiles_matching_subtype = int(st.n_tiles_matching_positive)
        tiles_matching_subtype_expected = int(st.n_tiles_matching_positive_expected)
        missing_subtype_tiles = tiles_matching_subtype_expected - tiles_matching_subtype

        negative_tiles_matching_subtype = int(st.n_tiles_matching_negative)
        negative_tiles_matching_subtype_expected = int(st.n_tiles_matching_negative_expected)
        missing_negative_subtype_tiles = negative_tiles_matching_subtype_expected - negative_tiles_matching_subtype

        total_missing_target_tiles = missing_subtype_tiles + missing_negative_subtype_tiles
        threshold_for_missing_tiles = int(st.n_tiles_matching_all_expected) - \
            (int(st.n_tiles_matching_all_expected) * MIN_TILES_THRESHOLD)

        if threshold_for_missing_tiles <= int(st.n_tiles_matching_all):
            if 3 <= total_missing_target_tiles and missing_subtype_tiles and missing_negative_subtype_tiles:
                error_status = FAIL_MESSAGE
                error_messages = ("{}: {} missing tiles detected for subtype: {}."
                                  " {} positive tiles missing, {} negative tiles missing."
                                  .format(
                                            INCONSISTENT_RESULTS_ERROR_3A,
                                            total_missing_target_tiles,
                                            st.subtype,
                                            missing_subtype_tiles,
                                            missing_negative_subtype_tiles
                                         )
                                  )
            possible_ds_subtypes = possible_subtypes_exist_in_df(st, df)
            if possible_ds_subtypes:
                error_status = FAIL_MESSAGE
                error_messages = ("{}: Subtype {} was found, but downstream subtype(s) {} tiles were missing."
                                  " Therefore due to the missing downstream tiles, there is a lack of confidence in the"
                                  "final subtype call."
                                  .format(
                                            INCONSISTENT_RESULTS_ERROR_3B,
                                            st.subtype,
                                            possible_ds_subtypes
                                         )
                                  )

    else:
        logging.debug("QC: Checking for inconsistent results not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for inconsistent results not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages


def check_intermediate_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    """ Check for intermediate subtypes within the result.
    Note:
            Check_intermediate_subtype will return a WARNING message if 95% of the scheme tiles are found, there are
            0 contradicting tiles (tiles with same position but there exists + and -), total subtype tiles < expected,
            and +/- tiles exist in the result.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Output:
            Tuple[Optional[str], Optional[str]]
            error_status: Contains any error statuses { Warning }
            error_messages: Contains the reasoning for any error status.
    """
    logging.debug("QC: Checking for intermediate subtypes.")
    error_status = None
    error_messages = None

    if st.are_subtypes_consistent:
        total_tiles = int(st.n_tiles_matching_all_expected)
        total_tiles_hits = int(st.n_tiles_matching_all)
        total_subtype_tiles = int(st.n_tiles_matching_subtype_expected)
        total_subtype_tiles_hits = int(st.n_tiles_matching_subtype)
        conflicting_tiles = get_conflicting_tiles(st, df)
        num_pos_tiles, num_neg_tiles = get_num_pos_neg_tiles(st, df)

        if (total_tiles - (total_tiles * MIN_TILES_THRESHOLD)) <= total_tiles_hits and not conflicting_tiles \
                and total_subtype_tiles_hits < total_subtype_tiles and num_pos_tiles and num_neg_tiles:

            error_status = WARNING_MESSAGE
            error_messages = "{}: Possible intermediate subtype. All scheme tiles were found, " \
                             "but a fraction were positive for the final subtype. " \
                             "Total subtype hits: {} | Total subtype expected: {}." \
                             .format(INTERMEDIATE_SUBTYPE_WARNING, total_subtype_tiles_hits, total_subtype_tiles)
    else:
        logging.debug("QC: Checking for intermediate subtypes not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for intermediate subtype not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages
