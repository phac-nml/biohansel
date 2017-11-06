from bio_hansel.const import MIXED_SUBTYPE_ERROR, INSUFFICIENT_NUM_TILES,\
    CONFIDENT_SUBTYPE_ERROR, FAIL_MESSAGE, PASS_MESSAGE, MIN_TILES_ERROR
from .subtype import Subtype


def perform_quality_check(st: Subtype):
    qc_message = ''

    if st is not None:
        qc_status_confidence, qc_message_confidence = check_is_confident_subtype(st)
        qc_status_min_tiles, qc_message_min_tiles = check_min_tiles_reached(st)

        if qc_status_confidence != qc_status_min_tiles and qc_status_confidence is FAIL_MESSAGE or qc_status_min_tiles \
                is FAIL_MESSAGE:
            qc_status = FAIL_MESSAGE
        else:
            qc_status = PASS_MESSAGE

        if 0 < len(qc_message_confidence):
            qc_message += qc_message_confidence

        if 0 < len(qc_message_min_tiles):
            if 0 < len(qc_message):
                qc_message += ' | '
            qc_message += qc_message_min_tiles

        st.qc_message = qc_message
        st.qc_status = qc_status


def check_is_confident_subtype(st: Subtype) -> str:
    # Value is default false, unless checks below pass which is where it's set to true.
    does_subtype_pass = FAIL_MESSAGE
    error_messages = ''
    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
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
    else:
        error_messages = CONFIDENT_SUBTYPE_ERROR
        does_subtype_pass = FAIL_MESSAGE

    return does_subtype_pass, error_messages


def check_min_tiles_reached(st: Subtype) -> str:
    does_subtype_pass = FAIL_MESSAGE
    error_messages = ''

    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
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
    else:
        does_subtype_pass = FAIL_MESSAGE
        error_messages = MIN_TILES_ERROR

    return does_subtype_pass, error_messages
