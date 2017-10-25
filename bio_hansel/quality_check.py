from bio_hansel.const import MIXED_SUBTYPE_ERROR, OK_SUBTYPE, INSUFFICIENT_NUM_TILES, OK_NUM_TILES
from .subtype import Subtype


def check_is_confident_subtype(st: Subtype) -> str:
    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
        if ';' in str(st.n_tiles_matching_all_expected):
            expected_tiles_matching = int(st.n_tiles_matching_all_expected.split(";")[0])
        else:
            expected_tiles_matching = int(st.n_tiles_matching_all_expected)

        if st.are_subtypes_consistent is False or (st.inconsistent_subtypes is not None and st.inconsistent_subtypes > 0):
            confident_subtype = MIXED_SUBTYPE_ERROR
        elif st.n_tiles_matching_all >= ((expected_tiles_matching * 0.01) +
                                         expected_tiles_matching):
            confident_subtype = MIXED_SUBTYPE_ERROR
        else:
            confident_subtype = OK_SUBTYPE
    else:
        confident_subtype = ""

    return confident_subtype


def check_min_tiles_reached(st: Subtype) -> str:
    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
        if ';' in str(st.n_tiles_matching_all_expected):
            expected_tiles_matching = int(st.n_tiles_matching_all_expected.split(";")[0])
        else:
            expected_tiles_matching = int(st.n_tiles_matching_all_expected)

        if st.n_tiles_matching_all <= expected_tiles_matching - (expected_tiles_matching * 0.05):
            min_tiles_reached = INSUFFICIENT_NUM_TILES
        else:
            min_tiles_reached = OK_NUM_TILES
    else:
        min_tiles_reached = ""

    return min_tiles_reached
