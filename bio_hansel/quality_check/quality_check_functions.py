from bio_hansel.quality_check.const import FAIL_MESSAGE, MIXED_SUBTYPE_ERROR, PASS_MESSAGE, INSUFFICIENT_NUM_TILES
from bio_hansel.subtype import Subtype


def does_subtype_result_exist(st):
    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
        return True
    else:
        return False


def check_is_confident_subtype(st: Subtype) -> str:
    # Value is default false, unless checks below pass which is where it's set to true.
    does_subtype_pass = FAIL_MESSAGE
    error_messages = ''
    # First we need to get the expected tiles that are matching, this is stored as a string.
    if ';' in str(st.n_tiles_matching_all_expected):
        expected_tiles_matching = int(st.n_tiles_matching_all_expected.split(";")[0])
    else:
        expected_tiles_matching = int(st.n_tiles_matching_all_expected)

        # Now we need to make sure that the subtypes are consistent or provide an explanation of what is wrong.
    if st.are_subtypes_consistent is False or (st.inconsistent_subtypes is not None and st.inconsistent_subtypes > 0):
        error_messages = MIXED_SUBTYPE_ERROR
    # Checking the max tiles threshold, to indicate whether there's a mixed subtype or not.
    elif st.n_tiles_matching_all >= ((expected_tiles_matching * 0.01) +
                                     expected_tiles_matching):
        error_messages = MIXED_SUBTYPE_ERROR
    else:
        does_subtype_pass = PASS_MESSAGE

    return does_subtype_pass, error_messages


def check_min_tiles_reached(st: Subtype) -> str:
    does_subtype_pass = FAIL_MESSAGE
    error_messages = ''

    if ';' in str(st.n_tiles_matching_all_expected):
        # Need to calculate the minimum tiles expected, which is stored as a string.
        expected_tiles_matching = int(st.n_tiles_matching_all_expected.split(";")[0])
    else:
        expected_tiles_matching = int(st.n_tiles_matching_all_expected)

    # Then we verify that the subtype has the correct number of tiles.
    if st.n_tiles_matching_all <= expected_tiles_matching - (expected_tiles_matching * 0.05):
        error_messages = INSUFFICIENT_NUM_TILES
    else:
        does_subtype_pass = PASS_MESSAGE

    return does_subtype_pass, error_messages

