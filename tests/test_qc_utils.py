# -*- coding: utf-8 -*-

import pandas as pd

from bio_hansel.qc.utils import get_mixed_subtype_kmer_counts


def test_get_mixed_subtype_kmer_counts():
    """
    Subtype 2.1 should have positive kmers (2, 2, 2.1, 2.1, 2.1) == 5 pos kmers
    Subtype 2.2 should have positive kmers (2, 2, 2.2) == 3 pos kmers
    Subtype 0 should have 1 positive kmer (0)
    """
    subtype_list = ['2.1', '2.2', '0']
    data = {'subtype': ['1', '2', '2', '2.1', '2.1', '2.1', '2.2', '2.3', '2.4', '0']}
    df = pd.DataFrame(data=data)

    st_pos_kmers = get_mixed_subtype_kmer_counts(df, subtype_list)

    assert (isinstance(st_pos_kmers, dict))
    assert(len(st_pos_kmers) == 3)
    assert(int(st_pos_kmers.get('2.1')) == 5)
    assert(int(st_pos_kmers.get('2.2')) == 3)
    assert(int(st_pos_kmers.get('0')) == 1)
