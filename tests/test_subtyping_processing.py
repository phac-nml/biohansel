# -*- coding: utf-8 -*-

import pytest
import pandas as pd

from bio_hansel.qc import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import absent_downstream_subtypes, sorted_subtype_ints, count_periods, empty_results, subtype_sets, missing_nested_subtypes
from bio_hansel.utils import find_inconsistent_subtypes
from bio_hansel.subtype_stats import SubtypeCounts
from bio_hansel.aho_corasick import expand_degenerate_bases


def test_count_periods():
    assert count_periods('1.1.1') == 2
    assert count_periods('1.1.10000') == 2
    assert count_periods('1.10000') == 1
    assert count_periods('10000') == 0


def test_absent_downstream_subtypes():
    assert absent_downstream_subtypes('1', pd.Series(['1.1', '1.2', '1.3', '1']), [
                                      '1.1', '1.2', '1', '1.3']) is None
    assert absent_downstream_subtypes('1', pd.Series(['1.1', '1.2', '1']), [
                                      '1.1', '1.2', '1', '1.3']) == ['1.3']
    assert absent_downstream_subtypes(
        '1', pd.Series(
            ['1']), [
            '1.1', '1.2', '1', '1.3']) == [
                '1.1', '1.2', '1.3']


def test_sorted_subtype_ints():
    assert sorted_subtype_ints(pd.Series([])) == []
    assert sorted_subtype_ints(pd.Series(['1', '1.1', '1.1.1', '1.1.1.99'])) == [
        [1], [1, 1], [1, 1, 1], [1, 1, 1, 99]]
    assert sorted_subtype_ints(
        pd.Series(
            [
                '1', '1.1', '1.1.1', '1.1.1.99', '1.1', '1.1.1'])) == [
        [1], [
            1, 1], [
            1, 1, 1], [
            1, 1, 1, 99]]


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
        'a']
    for good_value in good_values:
        assert SubtypeCounts._check_subtype('x', 'x', good_value) == good_value
    with pytest.raises(ValueError):
        for bad_value in bad_values:
            assert SubtypeCounts._check_subtype('x', 'x', bad_value) == ''


def test_subtype_sets():
    st_vals = ['2', '22', '222', '2222', '22222']
    pos_subtypes_set = {'2', '2.22.222.2222.22222'}
    assert subtype_sets(
        st_vals,
        pos_subtypes_set,
        primary_subtypes_set=set()) == {
        '2.22',
        '2.22.222',
        '2.22.222.2222'}


def test_expand_degenerate_bases():
    assert len(expand_degenerate_bases('NNNNN')) == 1024
    assert expand_degenerate_bases('DARTHVADR') == [
        'AAATAAAAA',
        'AAATAAAAG',
        'AAATAAAGA',
        'AAATAAAGG',
        'AAATAAATA',
        'AAATAAATG',
        'AAATACAAA',
        'AAATACAAG',
        'AAATACAGA',
        'AAATACAGG',
        'AAATACATA',
        'AAATACATG',
        'AAATAGAAA',
        'AAATAGAAG',
        'AAATAGAGA',
        'AAATAGAGG',
        'AAATAGATA',
        'AAATAGATG',
        'AAATCAAAA',
        'AAATCAAAG',
        'AAATCAAGA',
        'AAATCAAGG',
        'AAATCAATA',
        'AAATCAATG',
        'AAATCCAAA',
        'AAATCCAAG',
        'AAATCCAGA',
        'AAATCCAGG',
        'AAATCCATA',
        'AAATCCATG',
        'AAATCGAAA',
        'AAATCGAAG',
        'AAATCGAGA',
        'AAATCGAGG',
        'AAATCGATA',
        'AAATCGATG',
        'AAATTAAAA',
        'AAATTAAAG',
        'AAATTAAGA',
        'AAATTAAGG',
        'AAATTAATA',
        'AAATTAATG',
        'AAATTCAAA',
        'AAATTCAAG',
        'AAATTCAGA',
        'AAATTCAGG',
        'AAATTCATA',
        'AAATTCATG',
        'AAATTGAAA',
        'AAATTGAAG',
        'AAATTGAGA',
        'AAATTGAGG',
        'AAATTGATA',
        'AAATTGATG',
        'AAGTAAAAA',
        'AAGTAAAAG',
        'AAGTAAAGA',
        'AAGTAAAGG',
        'AAGTAAATA',
        'AAGTAAATG',
        'AAGTACAAA',
        'AAGTACAAG',
        'AAGTACAGA',
        'AAGTACAGG',
        'AAGTACATA',
        'AAGTACATG',
        'AAGTAGAAA',
        'AAGTAGAAG',
        'AAGTAGAGA',
        'AAGTAGAGG',
        'AAGTAGATA',
        'AAGTAGATG',
        'AAGTCAAAA',
        'AAGTCAAAG',
        'AAGTCAAGA',
        'AAGTCAAGG',
        'AAGTCAATA',
        'AAGTCAATG',
        'AAGTCCAAA',
        'AAGTCCAAG',
        'AAGTCCAGA',
        'AAGTCCAGG',
        'AAGTCCATA',
        'AAGTCCATG',
        'AAGTCGAAA',
        'AAGTCGAAG',
        'AAGTCGAGA',
        'AAGTCGAGG',
        'AAGTCGATA',
        'AAGTCGATG',
        'AAGTTAAAA',
        'AAGTTAAAG',
        'AAGTTAAGA',
        'AAGTTAAGG',
        'AAGTTAATA',
        'AAGTTAATG',
        'AAGTTCAAA',
        'AAGTTCAAG',
        'AAGTTCAGA',
        'AAGTTCAGG',
        'AAGTTCATA',
        'AAGTTCATG',
        'AAGTTGAAA',
        'AAGTTGAAG',
        'AAGTTGAGA',
        'AAGTTGAGG',
        'AAGTTGATA',
        'AAGTTGATG',
        'GAATAAAAA',
        'GAATAAAAG',
        'GAATAAAGA',
        'GAATAAAGG',
        'GAATAAATA',
        'GAATAAATG',
        'GAATACAAA',
        'GAATACAAG',
        'GAATACAGA',
        'GAATACAGG',
        'GAATACATA',
        'GAATACATG',
        'GAATAGAAA',
        'GAATAGAAG',
        'GAATAGAGA',
        'GAATAGAGG',
        'GAATAGATA',
        'GAATAGATG',
        'GAATCAAAA',
        'GAATCAAAG',
        'GAATCAAGA',
        'GAATCAAGG',
        'GAATCAATA',
        'GAATCAATG',
        'GAATCCAAA',
        'GAATCCAAG',
        'GAATCCAGA',
        'GAATCCAGG',
        'GAATCCATA',
        'GAATCCATG',
        'GAATCGAAA',
        'GAATCGAAG',
        'GAATCGAGA',
        'GAATCGAGG',
        'GAATCGATA',
        'GAATCGATG',
        'GAATTAAAA',
        'GAATTAAAG',
        'GAATTAAGA',
        'GAATTAAGG',
        'GAATTAATA',
        'GAATTAATG',
        'GAATTCAAA',
        'GAATTCAAG',
        'GAATTCAGA',
        'GAATTCAGG',
        'GAATTCATA',
        'GAATTCATG',
        'GAATTGAAA',
        'GAATTGAAG',
        'GAATTGAGA',
        'GAATTGAGG',
        'GAATTGATA',
        'GAATTGATG',
        'GAGTAAAAA',
        'GAGTAAAAG',
        'GAGTAAAGA',
        'GAGTAAAGG',
        'GAGTAAATA',
        'GAGTAAATG',
        'GAGTACAAA',
        'GAGTACAAG',
        'GAGTACAGA',
        'GAGTACAGG',
        'GAGTACATA',
        'GAGTACATG',
        'GAGTAGAAA',
        'GAGTAGAAG',
        'GAGTAGAGA',
        'GAGTAGAGG',
        'GAGTAGATA',
        'GAGTAGATG',
        'GAGTCAAAA',
        'GAGTCAAAG',
        'GAGTCAAGA',
        'GAGTCAAGG',
        'GAGTCAATA',
        'GAGTCAATG',
        'GAGTCCAAA',
        'GAGTCCAAG',
        'GAGTCCAGA',
        'GAGTCCAGG',
        'GAGTCCATA',
        'GAGTCCATG',
        'GAGTCGAAA',
        'GAGTCGAAG',
        'GAGTCGAGA',
        'GAGTCGAGG',
        'GAGTCGATA',
        'GAGTCGATG',
        'GAGTTAAAA',
        'GAGTTAAAG',
        'GAGTTAAGA',
        'GAGTTAAGG',
        'GAGTTAATA',
        'GAGTTAATG',
        'GAGTTCAAA',
        'GAGTTCAAG',
        'GAGTTCAGA',
        'GAGTTCAGG',
        'GAGTTCATA',
        'GAGTTCATG',
        'GAGTTGAAA',
        'GAGTTGAAG',
        'GAGTTGAGA',
        'GAGTTGAGG',
        'GAGTTGATA',
        'GAGTTGATG',
        'TAATAAAAA',
        'TAATAAAAG',
        'TAATAAAGA',
        'TAATAAAGG',
        'TAATAAATA',
        'TAATAAATG',
        'TAATACAAA',
        'TAATACAAG',
        'TAATACAGA',
        'TAATACAGG',
        'TAATACATA',
        'TAATACATG',
        'TAATAGAAA',
        'TAATAGAAG',
        'TAATAGAGA',
        'TAATAGAGG',
        'TAATAGATA',
        'TAATAGATG',
        'TAATCAAAA',
        'TAATCAAAG',
        'TAATCAAGA',
        'TAATCAAGG',
        'TAATCAATA',
        'TAATCAATG',
        'TAATCCAAA',
        'TAATCCAAG',
        'TAATCCAGA',
        'TAATCCAGG',
        'TAATCCATA',
        'TAATCCATG',
        'TAATCGAAA',
        'TAATCGAAG',
        'TAATCGAGA',
        'TAATCGAGG',
        'TAATCGATA',
        'TAATCGATG',
        'TAATTAAAA',
        'TAATTAAAG',
        'TAATTAAGA',
        'TAATTAAGG',
        'TAATTAATA',
        'TAATTAATG',
        'TAATTCAAA',
        'TAATTCAAG',
        'TAATTCAGA',
        'TAATTCAGG',
        'TAATTCATA',
        'TAATTCATG',
        'TAATTGAAA',
        'TAATTGAAG',
        'TAATTGAGA',
        'TAATTGAGG',
        'TAATTGATA',
        'TAATTGATG',
        'TAGTAAAAA',
        'TAGTAAAAG',
        'TAGTAAAGA',
        'TAGTAAAGG',
        'TAGTAAATA',
        'TAGTAAATG',
        'TAGTACAAA',
        'TAGTACAAG',
        'TAGTACAGA',
        'TAGTACAGG',
        'TAGTACATA',
        'TAGTACATG',
        'TAGTAGAAA',
        'TAGTAGAAG',
        'TAGTAGAGA',
        'TAGTAGAGG',
        'TAGTAGATA',
        'TAGTAGATG',
        'TAGTCAAAA',
        'TAGTCAAAG',
        'TAGTCAAGA',
        'TAGTCAAGG',
        'TAGTCAATA',
        'TAGTCAATG',
        'TAGTCCAAA',
        'TAGTCCAAG',
        'TAGTCCAGA',
        'TAGTCCAGG',
        'TAGTCCATA',
        'TAGTCCATG',
        'TAGTCGAAA',
        'TAGTCGAAG',
        'TAGTCGAGA',
        'TAGTCGAGG',
        'TAGTCGATA',
        'TAGTCGATG',
        'TAGTTAAAA',
        'TAGTTAAAG',
        'TAGTTAAGA',
        'TAGTTAAGG',
        'TAGTTAATA',
        'TAGTTAATG',
        'TAGTTCAAA',
        'TAGTTCAAG',
        'TAGTTCAGA',
        'TAGTTCAGG',
        'TAGTTCATA',
        'TAGTTCATG',
        'TAGTTGAAA',
        'TAGTTGAAG',
        'TAGTTGAGA',
        'TAGTTGAGG',
        'TAGTTGATA',
        'TAGTTGATG']
