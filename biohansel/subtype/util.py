from typing import Optional

import os

from biohansel.subtype.const import SCHEME_FASTAS
from biohansel.subtype.subtyping_params import SubtypingParams


def get_scheme_fasta(scheme: str) -> str:
    if scheme in SCHEME_FASTAS:
        scheme_fasta = SCHEME_FASTAS[scheme]['file']
    elif os.path.exists(scheme) and os.path.isfile(scheme):
        scheme_fasta = scheme
    else:
        raise FileNotFoundError('Could not find user-specified subtyping scheme fasta "%s"', scheme)
    return scheme_fasta


def get_scheme_params(scheme: str) -> Optional[SubtypingParams]:
    if scheme in SCHEME_FASTAS:
        return SCHEME_FASTAS[scheme]['subtyping_params']


def get_scheme_version(scheme: str) -> Optional[str]:
    if scheme in SCHEME_FASTAS:
        version = SCHEME_FASTAS[scheme]['version']  # type: str
        return version
    return None


def init_subtyping_params(low_coverage_threshold: float = None,
                          low_coverage_warning: float = None,
                          max_intermediate_tiles: float = None,
                          max_missing_tiles: float = None,
                          min_ambiguous_tiles: int = None,
                          scheme: str = None,
                          **kwargs) -> SubtypingParams:
    """Initialize subtyping parameters based on command-line arguments and scheme defaults

    Args:
        args: ArgumentParser.parse_args() output
        scheme: Scheme name e.g. "heidelberg"

    Returns:
        SubtypingParams with user-supplied values then scheme defaults then global defaults loaded
    """
    subtyping_params = get_scheme_params(scheme)
    if subtyping_params is None:
        subtyping_params = SubtypingParams()
    if low_coverage_threshold:
        subtyping_params.low_coverage_threshold = low_coverage_threshold
    if max_missing_tiles:
        subtyping_params.max_missing_tiles = max_missing_tiles
    if min_ambiguous_tiles:
        subtyping_params.min_ambiguous_tiles = min_ambiguous_tiles
    if max_intermediate_tiles:
        subtyping_params.max_intermediate_tiles = max_intermediate_tiles
    if low_coverage_warning:
        subtyping_params.low_coverage_warning = low_coverage_warning

    return subtyping_params
