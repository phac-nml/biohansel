from pandas import DataFrame

from bio_hansel.quality_check.const import FAIL_MESSAGE, MIXED_SUBTYPE_ERROR, INSUFFICIENT_NUM_TILES, \
    MAX_TILES_THRESHOLD, MIXED_SUBTYPE_WARNING, WARNING_MESSAGE, OVER_MAX_TILES, MIN_TILES_THRESHOLD
from bio_hansel.subtype import Subtype
from typing import Tuple, Optional


def check_missing_tiles(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    positive_tile_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])]
    positive_tile_totals = df[(df['subtype'] == st.subtype) & (df['is_pos_tile'])]

    negative_tile_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])]
    negative_tile_totals = df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False)]

    total_tiles = len(positive_tile_totals) + len(negative_tile_totals)

    if (len(positive_tile_hits) + len(negative_tile_hits)) < (total_tiles - total_tiles * MIN_TILES_THRESHOLD):

        error_status = FAIL_MESSAGE
        # TODO: Need to calculate genome coverage, display within this error message.
        error_messages = "Missing tiles were detected. Calculated "
    return True


def check_mixed_subtype(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = []

    positive_tile_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])]
    negative_tile_hits = df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])]

    if st.are_subtypes_consistent is False or len(st.inconsistent_subtypes) > 0:
        error_status = FAIL_MESSAGE
        error_messages.append("Inconsistent subtypes detected. Subtypes found: {}.".format(st.inconsistent_subtypes))
    if 0 < len(positive_tile_hits) and 0 < len(negative_tile_hits):
        error_status = FAIL_MESSAGE
        error_messages.append("Positive and negative tiles detected for the same subtype {}.".format(st.subtype))

    error_messages = ' | '.join(error_messages)
    return error_status, error_messages


def check_inconsistent_results(st: Subtype, df: DataFrame) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    positive_tile_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile']) & (df['is_kmer_freq_okay'])])
    positive_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'])])
    negative_tile_hits = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False) & (df['is_kmer_freq_okay'])])
    negative_tile_totals = len(df[(df['subtype'] == st.subtype) & (df['is_pos_tile'] == False)])

    total_tiles = positive_tile_totals + negative_tile_totals
    total_hits = positive_tile_hits + negative_tile_hits

    if (total_tiles - total_tiles * 0.05) <= total_hits <= (total_tiles + total_tiles * 0.05):
        total_missing_tiles = ((positive_tile_totals - positive_tile_hits) + (negative_tile_totals - negative_tile_hits))
        if 3 <= total_missing_tiles:
            error_status = FAIL_MESSAGE
            error_messages = ("Inconsistent Results: %d missing tiles detected for subtype: {}.".format
                                  (total_missing_tiles, st.subtype))

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
