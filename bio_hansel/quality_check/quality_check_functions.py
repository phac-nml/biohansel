from pandas import DataFrame

from bio_hansel.quality_check.const import FAIL_MESSAGE, MIXED_SUBTYPE_ERROR, INSUFFICIENT_NUM_TILES, \
    MAX_TILES_THRESHOLD, MIXED_SUBTYPE_WARNING, WARNING_MESSAGE, OVER_MAX_TILES, MIN_TILES_THRESHOLD
from bio_hansel.subtype import Subtype
from typing import Tuple, Optional

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
    return st.subtype is not None and len(st.subtype) > 0


def check_missing_tiles(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    positive_tile_target_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])])
    positive_tile_target_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'])])

    negative_tile_target_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])])
    negative_tile_target_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False)])

    total_target_tiles = positive_tile_target_totals + negative_tile_target_totals

    if (positive_tile_target_hits + negative_tile_target_hits) < (total_target_tiles - total_target_tiles * MIN_TILES_THRESHOLD):
        # Genevieve verified this is the expected value.
        tiles_with_hits = df[(df['is_kmer_freq_okay'] == True)]
        average_freq_coverage_depth = tiles_with_hits['freq'].sum() / len(tiles_with_hits)

        error_messages = "More than 5% missing tiles were detected." \
                         " Average calculated tile coverage = {}".format(str(average_freq_coverage_depth))
        error_status = FAIL_MESSAGE

    return error_status, error_messages


def check_mixed_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = []

    positive_tile_target_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])]
    negative_tile_target_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])]

    if st.are_subtypes_consistent is False or len(st.inconsistent_subtypes) > 0:
        error_status = FAIL_MESSAGE
        error_messages.append("Mixed subtypes detected. "
                              "Mixed subtypes found: {}.".format(st.inconsistent_subtypes))
    if 0 < len(positive_tile_target_hits) and 0 < len(negative_tile_target_hits):
        error_status = FAIL_MESSAGE
        error_messages.append("Positive and negative tiles detected for the same subtype {}.".format(st.subtype))

    error_messages = ' | '.join(error_messages)
    return error_status, error_messages


def check_inconsistent_results(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    positive_tile_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])])
    positive_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'])])
    negative_tile_hits = len(
        df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])])
    negative_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False)])

    if st.are_subtypes_consistent:
        total_missing_target_tiles = (
            (positive_tile_totals - positive_tile_hits) + (negative_tile_totals - negative_tile_hits))
        threshold_for_missing_tiles = (
            int(st.n_tiles_matching_all_expected) - (int(st.n_tiles_matching_all_expected) * 0.05))

        if 3 <= total_missing_target_tiles and \
                (0 < positive_tile_hits and 0 < negative_tile_hits) and threshold_for_missing_tiles <= int(st.n_tiles_matching_all):
            error_status = FAIL_MESSAGE
            error_messages = ("Inconsistent Results: %d missing tiles detected for subtype: {}.".format
                              (total_missing_target_tiles, st.subtype))

    return error_status, error_messages


'''
def check_intermediate_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    positive_tile_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])])
    positive_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'])])

    negative_tile_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])])
    negative_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False)])

    return True

    # check for 5% missing targets
    # check for missing positive matches for expected targets
'''
