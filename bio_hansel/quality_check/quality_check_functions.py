from bio_hansel.quality_check.const import FAIL_MESSAGE, MIXED_SUBTYPE_ERROR, INSUFFICIENT_NUM_TILES, \
    MAX_TILES_THRESHOLD, MIXED_SUBTYPE_WARNING, WARNING_MESSAGE, OVER_MAX_TILES, MIN_TILES_THRESHOLD
from bio_hansel.subtype import Subtype
from typing import Tuple, Optional


''' 
[does_subtype_result_exist]
Input: Subtype 
Output: Bool,
        True: If the subtype contains any matching tiles.
        False: If the subtype does not contain any matching tiles.
Desc: This method verifies that there is an expected matching subtype so we can proceed
with quality checking.
'''


def does_subtype_result_exist(st):
    if st.n_tiles_matching_all_expected is not None and len(st.n_tiles_matching_all_expected) > 0:
        return True
    else:
        return False


''' 
[check_is_consistent_subtype]
Input: Subtype 
Output: Tuple, containing two strings, 
        error_status: Contains the status of the error, otherwise returns None.
        error_messages: Contains the error messages, otherwise returns None.
Desc: This method will verify that the subtype given to it is consistent by checking
the consistent flag and if there are any inconsistent subtypes within the st object.
'''


def check_is_consistent_subtype(st: Subtype) -> Tuple[Optional[str], Optional[str]]:
    # Value is default false, unless checks below pass which is where it's set to true.
    error_status = None
    error_messages = None

    if st.are_subtypes_consistent is False and (st.inconsistent_subtypes is not None
                                                and len(st.inconsistent_subtypes) > 0):
        mixed_subtypes = '; '.join(st.inconsistent_subtypes)
        error_messages = MIXED_SUBTYPE_ERROR + ": {" + mixed_subtypes + "} detected in sample " \
                                "{" + st.sample + "}. A single subtype is expected. " \
                                "This result could be due to contamination resulting in a mixed sample."
        error_status = FAIL_MESSAGE

    return error_status, error_messages


''' 
[check_min_tiles_reached]
Input: Subtype 
Output: Tuple, containing two strings, 
        error_status: Contains the status of the error, otherwise returns None.
        error_messages: Contains the error messages, otherwise returns None.
Desc: Checks if the minimum tiles have been reached within the subtype.
The threshold for this is to see if the `tiles matching all` is at least within 5% of the
expected tiles matching value. This method will return a FAIL and the offending tiles explanations
if the minimum number of tiles were not reached. This method will return a WARNING and an explanation
that the method was not run, if mixed subtypes were detected. 
'''


def check_min_tiles_reached(st: Subtype) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    if ';' not in str(st.n_tiles_matching_all_expected):
        expected_tiles_matching = int(st.n_tiles_matching_all_expected)

        # Then we verify that the subtype has the correct number of tiles.
        if st.n_tiles_matching_all <= expected_tiles_matching - (expected_tiles_matching * MIN_TILES_THRESHOLD):
            error_messages = INSUFFICIENT_NUM_TILES + " : Observed: {"+str(st.n_tiles_matching_all)+"} " \
                                                        " Expected: {"+str(st.n_tiles_matching_all_expected)+"}"
            error_status = FAIL_MESSAGE
    else:
        error_messages = MIXED_SUBTYPE_WARNING + "{Min Tiles QC}"
        error_status = WARNING_MESSAGE

    return error_status, error_messages


''' 
[check_max_tiles_reached]
Input: Subtype 
Output: Tuple, containing two strings, 
        error_status: Contains the status of the error, otherwise returns None.
        error_messages: Contains the error messages, otherwise returns None.
Desc: Checks if the maximum number of tiles are exceeded ( 1% ), if this is true this method
will return a FAIL status, and the error message containing the offending values.
If this method sees that there was a mixed subtype, this method won't run it's checks
and simply return a warning as the checks were not finished.
'''


def check_max_tiles_reached(st: Subtype) -> Tuple[Optional[str], Optional[str]]:
    error_status = None
    error_messages = None

    if ';' not in str(st.n_tiles_matching_all_expected):
        expected_tiles_matching = int(st.n_tiles_matching_all_expected)

        if st.n_tiles_matching_all >= ((expected_tiles_matching * MAX_TILES_THRESHOLD) +
                                       expected_tiles_matching):
            error_messages = OVER_MAX_TILES + " : Observed: {"+str(st.n_tiles_matching_all)+"}" \
                                                " Expected: {"+str(st.n_tiles_matching_all_expected)+"}"
            error_status = FAIL_MESSAGE
    else:
        error_messages = MIXED_SUBTYPE_WARNING + "{Max Tiles QC}"
        error_status = WARNING_MESSAGE

    return error_status, error_messages
