from typing import Tuple

from bio_hansel.subtype import Subtype
from pandas import DataFrame

'''
[get_conflicting_tiles]
Input: Subtype, DataFrame
Output: List
Desc: The purpose of this method is to find positive and negative tiles for the same refposition in the dataframe.
The method will return a list with the conflicting tiles.
'''


def get_conflicting_tiles(st: Subtype, df: DataFrame) -> list:
    if st.subtype is not None:
        if 'is_kmer_freq_okay' in df:
            dfst = df[(df['subtype'] == str(st.subtype)) & (df['is_kmer_freq_okay'])]
        else:  # fasta files
            dfst = df[(df['subtype'] == str(st.subtype))]

        pos_tiles = dfst[dfst['is_pos_tile']]
        neg_tiles = dfst[dfst['is_pos_tile'] == False]
        pos_tile_values = '|'.join(pos_tiles['refposition'].values.tolist())
        conflicting_tiles = neg_tiles[neg_tiles['refposition'].str.contains(pos_tile_values)][
            'refposition'].values.tolist()

    return conflicting_tiles


'''
[get_num_pos_neg_tiles]
Input: Subtype, DataFrame
Output: Tuple[int,int]
Desc: The purpose of this method is to find the number of positive and negative tiles that exist for a subtype.
'''


def get_num_pos_neg_tiles(st: Subtype, df: DataFrame) -> Tuple[int, int]:
    num_pos_tiles = 0
    num_neg_tiles = 0

    if st.subtype is not None:
        dfst = df[(df['subtype'] == str(st.subtype))]
        num_pos_tiles = dfst[dfst['is_pos_tile']].shape[0]
        num_neg_tiles = dfst[dfst['is_pos_tile'] == False].shape[0]

    return num_pos_tiles, num_neg_tiles
