from typing import List, Callable, Tuple

from pandas import DataFrame

from ..quality_check.quality_check_functions import check_missing_tiles, does_subtype_result_exist, \
    check_mixed_subtype, check_intermediate_subtype, check_inconsistent_results
from ..quality_check.const import FAIL_MESSAGE, WARNING_MESSAGE
from ..subtype import Subtype
import logging


QC_FUNCS: List[Callable[[Subtype, DataFrame], Tuple[str, str]]] = \
[
    check_missing_tiles,
    check_mixed_subtype,
    check_inconsistent_results,
    check_intermediate_subtype
]


def perform_quality_check(st: Subtype, df: DataFrame):
    """ Driver method to call all quality checking functions and handle their responses.
    Note:
            This is the driver method for the quality check module. Every method within the QC_FUNCS list will be run
            with parameters ( SUBTYPE, DATAFRAME ). If a quality check module returns something other than None, then
            an Error, or Warning has occured.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            None, modifies the subtype with the result.
    """
    logging.debug("Performing Quality Checking")
    overall_qc_status = 'PASS'
    messages = []

    if does_subtype_result_exist(st) is False:
        logging.warning("QC: Quality checking not run, subtype result did not exist.")
        st.qc_status = 'FAIL'
        st.qc_message = 'FAIL: Subtype does not exist, quality checking was not run.'
        return None

    for func in QC_FUNCS:
        # Calls run_method to check that the qc function takes a Subtype, returns Tuple[Optional[str], Optional[str]]
        status, message = func(st, df)
        if status is None:
            # If quality check function passes, move on to the next.
            continue
        messages.append('{}: {}'.format(status, message))
        if status is FAIL_MESSAGE:
            overall_qc_status = FAIL_MESSAGE
        elif overall_qc_status != FAIL_MESSAGE and status == WARNING_MESSAGE:
            overall_qc_status = WARNING_MESSAGE

    st.qc_status = overall_qc_status
    st.qc_message = ' | '.join(messages)
    logging.debug("QC: Finished!")
