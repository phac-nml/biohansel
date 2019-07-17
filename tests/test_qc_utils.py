# -*- coding: utf-8 -*-

import pandas as pd


from bio_hansel.qc.utils import get_mixed_subtype_kmer_counts, component_subtypes, get_conflicting_kmers


fail_tsv = 'tests/data/qc/conflicting_subtypes/fail.tsv'
pass_tsv = 'tests/data/qc/conflicting_subtypes/pass.tsv'


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


def test_component_subtypes():
    assert list(component_subtypes('4.2.1.1')) == ['4', '4.2', '4.2.1', '4.2.1.1']
    assert list(component_subtypes('1')) == ['1']


def test_get_conflicting_kmers():
    df_pass = pd.read_csv(pass_tsv, sep='\t')
    df_pass_result = get_conflicting_kmers('4.1.2', df_pass, True)
    assert df_pass_result.shape[0] == 0, 'Must be no conflicting kmers'
    df_fail = pd.read_csv(fail_tsv, sep='\t')
    df_fail_result = get_conflicting_kmers('4.1.2', df_fail, True)
    df_fail_result.reset_index(inplace=True)
    assert df_fail_result.shape[0] == 1, 'Must be one conflicting kmer'
    assert df_fail_result.refposition[0] == 62657
    assert df_fail_result.subtype[0] == '4.1'
    assert df_fail_result.is_pos_kmer[0] == False
