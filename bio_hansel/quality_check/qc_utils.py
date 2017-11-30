from typing import Tuple

from pandas import DataFrame

from ..subtype import Subtype


def get_conflicting_tiles(st: Subtype, df: DataFrame) -> list:
    """ This method gets positive and negative tiles that both are present for a subtype.
    Note:
            The purpose of this method is to find positive and negative tiles for the same refposition in the DataFrame.
            The method will return a list with the conflicting tiles.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            DataFrame containing the conflicting positive and negative tiles.
    """
    if st.subtype:
        if 'is_kmer_freq_okay' in df:
            dfst = df[(df['subtype'] == str(st.subtype)) & (df['is_kmer_freq_okay'])]
        else:  # fasta files
            dfst = df[(df['subtype'] == str(st.subtype))]

    pos_tile_positions = dfst[dfst['is_pos_tile']]['refposition'].tolist()
    neg_tiles = dfst[~dfst['is_pos_tile']]
    conflicting_tiles = neg_tiles[neg_tiles['refposition'].isin(pos_tile_positions)]

    return conflicting_tiles


def get_num_pos_neg_tiles(st: Subtype, df: DataFrame) -> Tuple[int, int]:
    """ This method gets the number of positive and negative tiles.
    Note:
            The purpose of this method is to find the count of positive and negative tiles, and return them to the
            caller.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            Tuple[int,int] containing the count of positive and negative tiles.
    """
    num_pos_tiles = 0
    num_neg_tiles = 0

    if st.subtype:
        dfst = df[(df['subtype'] == str(st.subtype))]
        num_pos_tiles = dfst[dfst['is_pos_tile']].shape[0]
        num_neg_tiles = dfst[~dfst['is_pos_tile']].shape[0]

    return num_pos_tiles, num_neg_tiles


