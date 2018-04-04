from typing import Tuple, Optional, List, Any, Dict

from pandas import DataFrame

from ..subtype import Subtype


def get_conflicting_tiles(st: Subtype, df: DataFrame) -> Optional[DataFrame]:
    """ Get positive and negative tiles that both are present for a subtype.

    Find positive and negative tiles for the same refposition/target site in the results `df`.

    Args:
        st: Subtype results summary
        df: Detailed subtyping results

    Returns:
        DataFrame of conflicting positive and negative tiles
    """
    if st.is_fastq_input():
        dfst = df[(df['subtype'] == str(st.subtype)) & (df['is_kmer_freq_okay'])]
    else:  # fasta files
        dfst = df[(df['subtype'] == str(st.subtype))]

    pos_tile_positions = dfst[dfst['is_pos_tile']]['refposition']
    neg_tiles = dfst[~dfst['is_pos_tile']]
    return neg_tiles[neg_tiles['refposition'].isin(pos_tile_positions)]


def get_num_pos_neg_tiles(st: Subtype, df: DataFrame) -> Tuple[int, int]:
    """Get the number of positive and negative tiles for a given subtype

    Args:
        st: Subtype results summary
        df: Detailed subtyping results

    Returns:
        (number of positive tiles, number of negative tiles)
    """
    dfst = df[(df['subtype'] == str(st.subtype))]
    return dfst[dfst['is_pos_tile']].shape[0], dfst[~dfst['is_pos_tile']].shape[0]


def get_mixed_subtype_tile_counts(dfpos: DataFrame, subtype_list: List[Any]) -> Dict[str, int]:
    st_pos_tiles = {}

    for st in subtype_list:
        bs = dfpos.subtype == [st[:x] for x in dfpos.subtype.str.len()]
        st_pos_tiles[st] = bs.sum()

    return st_pos_tiles
