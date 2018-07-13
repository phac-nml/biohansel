# -*- coding: utf-8 -*-

import pandas as pd

from biohansel.subtype.qc.utils import get_mixed_subtype_tile_counts


def test_get_mixed_subtype_tile_counts():
    """
    Subtype 2.1 should have positive tiles (2, 2, 2.1, 2.1, 2.1) == 5 pos tiles
    Subtype 2.2 should have positive tiles (2, 2, 2.2) == 3 pos tiles
    Subtype 0 should have 1 positive tile (0)
    """
    subtype_list = ['2.1', '2.2', '0']
    data = {'subtype': ['1', '2', '2', '2.1', '2.1', '2.1', '2.2', '2.3', '2.4', '0']}
    df = pd.DataFrame(data=data)

    st_pos_tiles = get_mixed_subtype_tile_counts(df, subtype_list)

    assert (isinstance(st_pos_tiles, dict))
    assert(len(st_pos_tiles) == 3)
    assert(int(st_pos_tiles.get('2.1')) == 5)
    assert(int(st_pos_tiles.get('2.2')) == 3)
    assert(int(st_pos_tiles.get('0')) == 1)
