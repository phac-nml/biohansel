from typing import Tuple, Optional, List, Any, Dict, Iterable

from pandas import DataFrame

from ..subtype import Subtype


def component_subtypes(subtype: str) -> Iterable[str]:
    """Generate component subtypes from a subtype.

    Args:
        subtype: Subtype string, e.g. "4.2.1.1"
    Yields:
        Component subtypes (e.g. for subtype "4.2.1.1", will yield
        ['4', '4.2', '4.2.1', '4.2.1.1'])
    """
    split_subtype = subtype.split('.')
    for i, x in enumerate(split_subtype):
        yield '.'.join(split_subtype[:i + 1])


def get_conflicting_kmers(subtype: str, df: DataFrame, is_fastq_input: bool = True) -> Optional[DataFrame]:
    """ Get positive and negative kmers that both are present for a subtype.

    Find positive and negative kmers for the same refposition/target site in the results `df`.

    Args:
        subtype: Subtype results summary
        df: Detailed subtyping results
        is_fastq_input: FASTQ input?

    Returns:
        DataFrame of conflicting positive and negative kmers
    """
    dfst = df[(df['subtype'].isin(list(component_subtypes(subtype))))]
    if is_fastq_input:
        dfst = dfst[dfst['is_kmer_freq_okay']]

    pos_kmer_positions = dfst[dfst['is_pos_kmer']]['refposition']
    neg_kmers = dfst[~dfst['is_pos_kmer']]
    return neg_kmers[neg_kmers['refposition'].isin(pos_kmer_positions)]


def get_num_pos_neg_kmers(st: Subtype, df: DataFrame) -> Tuple[int, int]:
    """Get the number of positive and negative kmers for a given subtype

    Args:
        st: Subtype results summary
        df: Detailed subtyping results

    Returns:
        (number of positive kmers, number of negative kmers)
    """
    dfst = df[(df['subtype'] == str(st.subtype))]
    return dfst[dfst['is_pos_kmer']].shape[0], dfst[~dfst['is_pos_kmer']].shape[0]


def get_mixed_subtype_kmer_counts(dfpos: DataFrame, subtype_list: List[Any]) -> Dict[str, int]:
    st_pos_kmers = {}

    for st in subtype_list:
        bs = dfpos.subtype == [st[:x] for x in dfpos.subtype.str.len()]
        st_pos_kmers[st] = bs.sum()

    return st_pos_kmers
