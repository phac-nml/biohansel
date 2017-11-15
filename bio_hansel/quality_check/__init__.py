from typing import List, Callable, Tuple, Optional

from pandas import DataFrame

from bio_hansel.quality_check.qc_utils import get_hit_and_miss_tiles
from bio_hansel.quality_check.quality_check_functions import check_missing_tiles, does_subtype_result_exist, \
    check_mixed_subtype, check_intermediate_subtype, check_inconsistent_results
from bio_hansel.quality_check.const import FAIL_MESSAGE, WARNING_MESSAGE
from bio_hansel.subtype import Subtype


QC_FUNCS: List[Callable[[Subtype, DataFrame], Tuple[str, str]]] = \
[
    check_missing_tiles,
    check_mixed_subtype,
    check_inconsistent_results,
    check_intermediate_subtype
]


'''
[perform_quality_check]
Input: Subtype, DataFrame
Output: None, modifies the subtype with the result.
Desc: This is the driver method for the quality check module. Every method within the QC_FUNCS list will be run 
with parameters ( SUBTYPE, DATAFRAME ). If a quality check module returns something other than None, then an Error,
or Warning has occured. 
'''


def perform_quality_check(st: Subtype, df: DataFrame):
    overall_qc_status = 'PASS'
    messages = []

    if does_subtype_result_exist(st) is False:
        st.qc_status = 'FAIL'
        st.qc_message = 'FAIL: Subtype does not exist, quality checking was not run.'
        return None

    for func in QC_FUNCS:
        # Calls run_method to check that the qc function takes a Subtype, returns Tuple[Optional[str], Optional[str]]
        status, message = run_method(func, st, df)
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


'''
[run_method]
Input: Function which takes [Subtype, DataFrame] outputs str tuple, Subtype, DataFrame
Output: Results of function Tuple[Optional[str], Optional[str]]
Desc: Helper method used to ensure function being called is of the right type. 
'''


def run_method(func: Callable[[Subtype, DataFrame], Tuple[str, str]], st: Subtype, df: DataFrame)\
        -> Tuple[Optional[str], Optional[str]]:
    return func(st, df)
