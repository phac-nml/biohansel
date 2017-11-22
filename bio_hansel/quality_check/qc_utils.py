from typing import Tuple
from ..subtype import Subtype
from pandas import DataFrame, to_numeric, Series


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
    dfst = df.copy()
    dfst['refposition'] = Series(dfst['refposition']).str.replace('negative', '')
    dfst['refposition'] = to_numeric(dfst['refposition'], downcast='unsigned', errors='coerce')

    if st.subtype:
        if 'is_kmer_freq_okay' in df:
            dfst = dfst[(dfst['subtype'] == str(st.subtype)) & (dfst['is_kmer_freq_okay'])]
        else:  # fasta files
            dfst = dfst[(dfst['subtype'] == str(st.subtype))]

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


def possible_subtypes_exist_in_df(st: Subtype, df: DataFrame) -> list:
    """ This method checks if the downstream subtypes' tiles are present within the DataFrame
    Note:
            The purpose of this method is to check if the downstream subtypes' tiles exist within the result.
            If they're not present then we know that the result may or not be confident.

    Args:
            :param st: Subtyping results.
            :param df: DataFrame containing subtyping results.

    Returns:
            list containing the non present subtypes.
    """
    non_present_subtypes = []
    possible_subtypes = st.possible_downstream_subtypes

    if possible_subtypes:
        for subtype in possible_subtypes:
            if subtype not in df['subtype']:
                non_present_subtypes.append(subtype)

    return non_present_subtypes
