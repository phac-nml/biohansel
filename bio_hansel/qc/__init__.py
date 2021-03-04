# -*- coding: utf-8 -*-

from typing import List, Callable, Tuple

from pandas import DataFrame

from ..qc.checks import \
    is_missing_kmers, \
    is_mixed_subtype, \
    is_maybe_intermediate_subtype, \
    is_missing_too_many_target_sites, \
    is_missing_downstream_targets, \
    is_missing_hierarchical_kmers, \
    is_overall_coverage_low
from ..qc.const import QC
from ..subtype import Subtype
from ..subtyping_params import SubtypingParams

CHECKS: List[Callable[[Subtype, DataFrame, SubtypingParams], Tuple[str, str]]] = [
    is_missing_kmers,
    is_mixed_subtype,
    is_missing_too_many_target_sites,
    is_missing_downstream_targets,
    is_missing_hierarchical_kmers,
    is_maybe_intermediate_subtype,
    is_overall_coverage_low
]


def perform_quality_check(st: Subtype, df: DataFrame, subtyping_params: SubtypingParams) -> Tuple[str, str]:
    """Perform QC of subtyping results

    Return immediate fail if subtype result is missing or if there are no detailed subtyping results.

    Args:
        st: Subtyping results.
        df: DataFrame containing subtyping results.
        subtyping_params: Subtyping/QC parameters

    Returns:
        (QC status, QC messages)
    """
    if st.subtype is None or len(st.subtype) == 0 \
            or df is None or df.shape[0] == 0:
        return QC.FAIL, QC.NO_SUBTYPE_RESULT

    overall_qc_status = QC.PASS
    messages = []
    for func in CHECKS:
        status, message = func(st, df, subtyping_params)
        # If quality check function passes, move on to the next.
        if status is None:
            continue
        messages.append('{}: {}'.format(status, message))
        if overall_qc_status == QC.FAIL:
            continue
        if status == QC.FAIL:
            overall_qc_status = QC.FAIL
            continue
        if status == QC.WARNING:
            overall_qc_status = QC.WARNING

    return overall_qc_status, ' | '.join(messages)
