from typing import Tuple, Optional

from pandas import DataFrame


def get_hit_and_miss_tiles(subtype: str, df: DataFrame) -> Tuple[Optional[int], Optional[int]]:
    hit = len(df[(df['subtype'] == subtype) & (df['is_pos_tile'])])
    miss = len(df[(df['subtype'] == subtype) & (df['is_pos_tile'] is False)])

    return hit, miss


def get_df_with_subtype(subtype: str, df: DataFrame) -> DataFrame:
    return df.loc[df['subtype'] == subtype]


def get_df_containing_subtype(subtype: str, df: DataFrame) -> DataFrame:
    return df.loc[df['subtype'].str.contains(subtype)]
