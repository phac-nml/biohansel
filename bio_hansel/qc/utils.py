from typing import Tuple, Optional, List, Any, Dict

from pandas import DataFrame

from ..subtype import Subtype


def get_conflicting_kmers(st: Subtype, df: DataFrame) -> Optional[DataFrame]:
    """ Get positive and negative kmers that both are present for a subtype.

    Find positive and negative kmers for the same refposition/target site in the results `df`.

    Args:
        st: Subtype results summary
        df: Detailed subtyping results

    Returns:
        DataFrame of conflicting positive and negative kmers
    """
    if st.is_fastq_input():
        dfst = df[(df['subtype'] == str(st.subtype)) & (df['is_kmer_freq_okay'])]
    else:  # fasta files
        dfst = df[(df['subtype'] == str(st.subtype))]

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
