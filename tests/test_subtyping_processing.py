# -*- coding: utf-8 -*-

import pytest
import pandas as pd

from bio_hansel.qc import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import absent_downstream_subtypes, sorted_subtype_ints, count_periods, empty_results
from bio_hansel.utils import find_inconsistent_subtypes
from bio_hansel.subtype_stats import SubtypeCounts


def test_count_periods():
    assert count_periods('1.1.1') == 2
    assert count_periods('1.1.10000') == 2
    assert count_periods('1.10000') == 1
    assert count_periods('10000') == 0


def test_absent_downstream_subtypes():
    assert absent_downstream_subtypes('1', pd.Series(['1.1', '1.2', '1.3', '1']), ['1.1', '1.2', '1', '1.3']) is None
    assert absent_downstream_subtypes('1', pd.Series(['1.1', '1.2', '1']), ['1.1', '1.2', '1', '1.3']) == ['1.3']
    assert absent_downstream_subtypes('1', pd.Series(['1']), ['1.1', '1.2', '1', '1.3']) == ['1.1', '1.2', '1.3']


def test_sorted_subtype_ints():
    assert sorted_subtype_ints(pd.Series([])) == []
    assert sorted_subtype_ints(pd.Series(['1', '1.1', '1.1.1', '1.1.1.99'])) == [[1], [1, 1], [1, 1, 1], [1, 1, 1, 99]]
    assert sorted_subtype_ints(pd.Series(['1', '1.1', '1.1.1', '1.1.1.99', '1.1', '1.1.1'])) == [[1], [1, 1], [1, 1, 1],
                                                                                                 [1, 1, 1, 99]]


def test_empty_results():
    st = Subtype(sample='test',
                 file_path='tests/data/Retro1000data/10-1358.fastq',
                 scheme='enteritidis',
                 scheme_version='0.8.0',
                 subtype=None,
                 non_present_subtypes=None,
                 all_subtypes=None,
                 qc_status=QC.FAIL,
                 qc_message=QC.NO_TARGETS_FOUND)
    df_empty = empty_results(st)
    df_expected_empty = pd.DataFrame({0: dict(sample='test',
                                              file_path='tests/data/Retro1000data/10-1358.fastq',
                                              subtype=None,
                                              refposition=None,
                                              is_pos_tile=None,
                                              scheme='enteritidis',
                                              scheme_version='0.8.0',
                                              qc_status=QC.FAIL,
                                              qc_message=QC.NO_TARGETS_FOUND)}).transpose()
    assert ((df_empty == df_expected_empty) | (df_empty.isnull() == df_expected_empty.isnull())).values.all(), \
        f'Empty result DataFrame should equal df_expected_empty: {df_expected_empty}'


def test_find_inconsistent_subtypes():
    consistent_subtypes = sorted_subtype_ints(pd.Series('''
1
1.1
1.1.1
1.1.1.1
    '''.strip().split('\n')))

    assert find_inconsistent_subtypes(consistent_subtypes) == [], \
        'Expecting all subtypes to be consistent with each other'

    inconsistent_subtypes = sorted_subtype_ints(pd.Series('''
1
1.1
1.1.1
1.1.1.1
1.1.1.2
1.1.1.3
    '''.strip().split('\n')))

    exp_incon_subtypes = '''
1.1.1.1
1.1.1.2
1.1.1.3
    '''.strip().split('\n')

    assert find_inconsistent_subtypes(inconsistent_subtypes) == exp_incon_subtypes, \
        f'Expecting subtypes {exp_incon_subtypes} to be inconsistent with each other'

    subtypes_list = '''
1
1.1
1.1.1
1.1.1.1
1.1.1.2
1.1.1.3
1.1.2
2
    '''.strip().split('\n')

    inconsistent_subtypes = sorted_subtype_ints(pd.Series(subtypes_list))
    assert set(find_inconsistent_subtypes(inconsistent_subtypes)) == set(subtypes_list), \
        f'All subtypes should be inconsistent with each other in {subtypes_list}'

def test_subtype_regex():
    good_values = ['1.1.1.1', '10', '77.10.1.9', '17.1.1.1.1.12.4']
    bad_values = ['1..', '1..1', '1.1..1.1', '1....', '100.', '', ' ', 'a1.1.1', '1.11.1a', 'a']
    for good_value in good_values:
        assert SubtypeCounts._check_subtype('x', 'x', good_value) == good_value
    with pytest.raises(ValueError):
        for bad_value in bad_values:
            assert SubtypeCounts._check_subtype('x', 'x', bad_value) == ''
        