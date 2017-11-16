from pandas import DataFrame
from bio_hansel.quality_check.const import FAIL_MESSAGE, WARNING_MESSAGE, MIN_TILES_THRESHOLD
from bio_hansel.quality_check.qc_utils import get_conflicting_tiles, get_num_pos_neg_tiles
from bio_hansel.subtype import Subtype
from typing import Tuple, Optional
import logging

''' 
[does_subtype_result_exist]
Input: Subtype 
Output: Bool,
        True: If the subtype exists
        False: If the subtype does not exist
Desc: This method verifies that there is an expected matching subtype so we can proceed
with quality checking.
'''


def does_subtype_result_exist(st) -> bool:
    logging.debug("QC: Checking if subtype result exists.")
    return st.subtype is not None and len(st.subtype) > 0


'''
[check_missing_tiles]
Input: Subtype, DataFrame
Output: Tuple[Optional[str], Optional[str]]
        error_status: Contains any error statuses { Warning, Fail }
        error_messages: Contains the reasoning for any error status.
Desc: check_missing_tiles will check if more than 5% of the scheme's tiles are missing.
If they are, we calculate the average frequency coverage depth for the tiles that are present, and
provide an adequate error message based on the coverage.
'''


def check_missing_tiles(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
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
                    error_messages = "More than 5% missing tiles were detected." \
                                     " Low coverage detected, possibly need more whole genome sequencing data." \
                                     " Average calculated tile coverage = {}".format(str(average_freq_coverage_depth))
                else:
                    error_messages = "More than 5% missing tiles were detected. " \
                                     "Adequate coverage detected, this may be the wrong serovar/species for scheme: {}"\
                            " Average calculated tile coverage = {}".format(str(st.scheme), average_freq_coverage_depth)
            else:
                error_messages = "More than 5% missing tiles were detected." \
                                 " Average calculated tile coverage = {}".format(str(average_freq_coverage_depth))

            error_status = FAIL_MESSAGE
    else:
        logging.debug("QC: Checking for missing tiles not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for missing tiles not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages


'''
[check_mixed_subtype]
Input: Subtype, DataFrame
Output: Tuple[Optional[str], Optional[str]]
        error_status: Contains any error statuses { Warning, Fail }
        error_messages: Contains the reasoning for any error status.
Desc: check_mixed_subtype will check the results of the typing for:
1) + and - tiles that are present for the target site above minimum coverage threshold.
2) If the subtyping result came back with an inconsistent call.
Then the method will provide an adequate error message to the user. 
'''


def check_mixed_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    logging.debug("QC: Checking for mixed subtypes.")
    error_status = None
    error_messages = None

    # First check the obvious, if the st is consistent, of course it's a mixed subtype.
    if st.are_subtypes_consistent is False or st.inconsistent_subtypes is not None:
        error_status = FAIL_MESSAGE
        error_messages = "Mixed subtypes detected. "\
                         "Mixed subtypes found: {}.".format(st.inconsistent_subtypes)
    else:
        # Find the conflicting tiles
        conflicting_tiles = get_conflicting_tiles(st, df)
        if 0 < len(conflicting_tiles):
            error_status = FAIL_MESSAGE
            error_messages = "Mixed subtype detected." \
                             " Positive and negative tiles detected for the same target site" \
                             " {} for subtype {}.".format(conflicting_tiles, st.subtype)

    return error_status, error_messages


'''
[check_inconsistent_results]
Input: Subtype, DataFrame
Output: Tuple[Optional[str], Optional[str]]
        error_status: Contains any error statuses { Warning, Fail }
        error_messages: Contains the reasoning for any error status.
Desc: An inconsistent results error occurs when 3 or more tiles are missing (+ / -) for the subtype call,
and there must be <= 5% tiles missing from the subtyping scheme.
This method returns the error_status == FAIL if inconsistent results were found, and the error_messages corresponding
to the error.
'''


def check_inconsistent_results(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    logging.debug("QC: Checking for inconsistent results.")
    error_status = None
    error_messages = None

    # We do this check as inconsistent subtypes fall within mixed subtypes. This checks for inconsistent results.
    if st.are_subtypes_consistent:

        positive_matching_hits = int(st.n_tiles_matching_positive)
        positive_matching_total = int(st.n_tiles_matching_positive_expected)
        positive_tile_misses = positive_matching_total - positive_matching_hits

        negative_matching_hits = int(st.n_tiles_matching_negative)
        negative_matching_total = int(st.n_tiles_matching_negative_expected)
        negative_tile_misses = negative_matching_total - negative_matching_hits

        total_missing_target_tiles = positive_tile_misses + negative_tile_misses
        threshold_for_missing_tiles = int(st.n_tiles_matching_all_expected) - (int(st.n_tiles_matching_all_expected) * MIN_TILES_THRESHOLD)

        if 3 <= total_missing_target_tiles and (0 < positive_tile_misses and 0 < negative_tile_misses) \
                and threshold_for_missing_tiles <= int(st.n_tiles_matching_all):
            error_status = FAIL_MESSAGE
            error_messages = ("Inconsistent Results: {} missing tiles detected for subtype: {}.".format
                              (total_missing_target_tiles, st.subtype))
    else:
        logging.debug("QC: Checking for inconsistent results not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for inconsistent results not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages


'''
[check_intermediate_subtype]
Input: Subtype
Output: Tuple[Optional[str], Optional[str]]
        error_status: Contains any error statuses { Warning }
        error_messages: Contains the reasoning for any error status.
Desc: check_intermediate_subtype will return a WARNING message if 95% of the scheme tiles are found, there are 
0 contradicting tiles (tiles with same position but there exists + and -), total subtype tiles < expected, and
+/- tiles exist in the result.
'''


def check_intermediate_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
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

        if (total_tiles - (total_tiles * MIN_TILES_THRESHOLD)) <= total_tiles_hits and len(conflicting_tiles) == 0 \
                and total_subtype_tiles_hits < total_subtype_tiles and (0 < num_pos_tiles and 0 < num_neg_tiles):

            error_status = WARNING_MESSAGE
            error_messages = "Possible intermediate subtype. All scheme tiles were found, " \
                             "but a fraction were positive for the final subtype. " \
                             "Total subtype hits: {} | Total subtype expected: {}." \
                             .format(total_subtype_tiles_hits, total_subtype_tiles)
    else:
        logging.debug("QC: Checking for intermediate subtypes not run, inconsistent subtype detected.")
        error_messages = "Subtype is inconsistent, quality checking for intermediate subtype not run."
        error_status = FAIL_MESSAGE

    return error_status, error_messages
