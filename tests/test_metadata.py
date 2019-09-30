# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

from bio_hansel.metadata import read_metadata_table, merge_results_with_metadata
from bio_hansel.utils import df_field_fillna


def test_read_metadata():
    df_results = pd.read_table('tests/data/subtyping-results.tsv')
    df_results = df_field_fillna(df_results)
    df_md = read_metadata_table('tests/data/subtype-metadata.tsv')
    df_merged = merge_results_with_metadata(df_results, df_md)
    assert np.all(df_md.columns.isin(df_merged.columns))
    for idx, row in df_merged.iterrows():
        # for metadata columns a,b,c, the merged metadata value should be the column name + subtype
        # e.g. for subtype 1.1, a=a1.1; b=b1.1; c=c1.1
        # if there was no subtype result, then the metadata field value should be NaN/None
        for column in ['a', 'b', 'c']:
            if row['subtype'] != '#N/A':
                assert row[column] == column + row['subtype']
            else:
                assert np.isnan(row[column])


def test_read_metadata_with_na_subtype():
    df_results = pd.read_table('tests/data/subtyping-results.tsv')
    df_results = df_field_fillna(df_results)
    df_md = read_metadata_table('tests/data/subtype-metadata-with-na.tsv')
    df_merged = merge_results_with_metadata(df_results, df_md)
    assert np.all(df_md.columns.isin(df_merged.columns))
    for idx, row in df_merged.iterrows():
        # for metadata columns a,b,c, the merged metadata value should be the column name + subtype
        # e.g. for subtype 1.1, a=a1.1; b=b1.1; c=c1.1
        # if there was no subtype result, then the metadata field value should be NaN/None
        for column in ['a', 'b', 'c']:
            if row['subtype'] != '#N/A':
                assert row[column] == column + row['subtype']
            else:
                assert row[column] == column, f'row={row}'
