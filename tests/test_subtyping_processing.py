# -*- coding: utf-8 -*-

import pandas as pd
import pytest

from bio_hansel.qc import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtype_stats import SubtypeCounts
from bio_hansel.subtyper import absent_downstream_subtypes, sorted_subtype_ints, empty_results, \
    get_missing_internal_subtypes
from bio_hansel.utils import find_inconsistent_subtypes, expand_degenerate_bases


def test_absent_downstream_subtypes():
    assert absent_downstream_subtypes(subtype='1',
                                      subtypes=pd.Series(['1.1', '1.2', '1.3', '1']),
                                      scheme_subtypes=['1.1', '1.2', '1', '1.3']) is None
    assert absent_downstream_subtypes(subtype='1',
                                      subtypes=pd.Series(['1.1', '1.2', '1']),
                                      scheme_subtypes=['1.1', '1.2', '1', '1.3']) == ['1.3']
    assert absent_downstream_subtypes(subtype='1',
                                      subtypes=pd.Series(['1']),
                                      scheme_subtypes=['1.1', '1.2', '1', '1.3']) == ['1.1', '1.2', '1.3']


def test_sorted_subtype_ints():
    assert sorted_subtype_ints(pd.Series([], dtype=object)) == []
    exp_subtype_ints = [
        [1],
        [1, 1],
        [1, 1, 1],
        [1, 1, 1, 99]
    ]
    assert sorted_subtype_ints(pd.Series(['1', '1.1', '1.1.1', '1.1.1.99'])) == exp_subtype_ints
    series = pd.Series(['1', '1.1', '1.1.1', '1.1.1.99', '1.1', '1.1.1'])
    assert sorted_subtype_ints(series) == exp_subtype_ints


def test_empty_results():
    st = Subtype(sample='test',
                 file_path='tests/data/Retro1000data/10-1358.fastq',
                 scheme='enteritidis',
                 scheme_version='1.0.5',
                 subtype=None,
                 non_present_subtypes=None,
                 all_subtypes=None,
                 qc_status=QC.FAIL,
                 qc_message=QC.NO_TARGETS_FOUND)
    df_empty = empty_results(st)
    df_expected_empty = pd.DataFrame(
        {
            0: dict(
                sample='test',
                file_path='tests/data/Retro1000data/10-1358.fastq',
                subtype=None,
                refposition=None,
                is_pos_kmer=None,
                scheme='enteritidis',
                scheme_version='1.0.5',
                qc_status=QC.FAIL,
                qc_message=QC.NO_TARGETS_FOUND)}).transpose()
    assert ((df_empty == df_expected_empty) | (df_empty.isnull() == df_expected_empty.isnull())).values.all(), \
        f'Empty result DataFrame should equal df_expected_empty: {df_expected_empty}'


def test_find_inconsistent_subtypes():
    subtype_list = ['1',
                    '1.1',
                    '1.1.1',
                    '1.1.1.1', ]
    consistent_subtypes = sorted_subtype_ints(pd.Series(subtype_list))

    assert find_inconsistent_subtypes(consistent_subtypes) == [], \
        'Expecting all subtypes to be consistent with each other'

    subtype_list = ['1',
                    '1.1',
                    '1.1.1',
                    '1.1.1.1',
                    '1.1.1.2',
                    '1.1.1.3', ]
    inconsistent_subtypes = sorted_subtype_ints(pd.Series(subtype_list))

    exp_incon_subtypes = ['1.1.1.1',
                          '1.1.1.2',
                          '1.1.1.3', ]
    assert find_inconsistent_subtypes(inconsistent_subtypes) == exp_incon_subtypes, \
        f'Expecting subtypes {exp_incon_subtypes} to be inconsistent with each other'

    subtypes_list = ['1',
                     '1.1',
                     '1.1.1',
                     '1.1.1.1',
                     '1.1.1.2',
                     '1.1.1.3',
                     '1.1.2',
                     '2', ]

    inconsistent_subtypes = sorted_subtype_ints(pd.Series(subtypes_list))
    assert set(find_inconsistent_subtypes(inconsistent_subtypes)) == set(subtypes_list), \
        f'All subtypes should be inconsistent with each other in {subtypes_list}'


def test_subtype_regex():
    good_values = ['1.1.1.1', '10', '77.10.1.9', '17.1.1.1.1.12.4', ]
    for good_value in good_values:
        assert SubtypeCounts._check_subtype(None, None, good_value) == good_value

    bad_values = [
        '1..',
        '1..1',
        '1.1..1.1',
        '1....',
        '100.',
        '',
        ' ',
        'a1.1.1',
        '1.11.1a',
        'a',
        'not.a.valid.subtype',
        'B.1.1.7'
    ]
    for bad_value in bad_values:
        with pytest.raises(ValueError):
            assert SubtypeCounts._check_subtype(None, None, bad_value) == ''


def test_get_missing_internal_subtypes():
    st_vals = ['1', '1', '1', '1']
    pos_subtypes_set = {
        '1',
        '1.1',
        '1.1.1',
        '1.1.1.1'
    }
    exp_missing_internal_subtypes = set()
    assert get_missing_internal_subtypes(st_vals, pos_subtypes_set) == exp_missing_internal_subtypes
    st_vals = ['2', '22', '222', '2222', '22222']
    pos_subtypes_set = {'2', '2.22.222.2222.22222'}
    exp_missing_internal_subtypes = {
        '2.22',
        '2.22.222',
        '2.22.222.2222'
    }
    assert get_missing_internal_subtypes(st_vals, pos_subtypes_set) == exp_missing_internal_subtypes


def test_expand_degenerate_bases():
    assert len(expand_degenerate_bases('NNNNN')) == 1024
    with open('tests/data/expand_degenerate_bases_DARTHVADR.txt') as f:
        assert expand_degenerate_bases('DARTHVADR') == f.read().split('\n')
