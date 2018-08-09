# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from biohansel.subtype.metadata import read_metadata_table, merge_metadata_with_summary_results


def test_read_metadata():
    df_results = pd.read_table('tests/data/subtype/subtyping-results.tsv')
    df_md = read_metadata_table('tests/data/subtype/subtype-metadata.tsv')
    df_merged = merge_metadata_with_summary_results(df_results, df_md)
    assert np.all(df_md.columns.isin(df_merged.columns))
    for _, row in df_merged.iterrows():
        # for metadata columns a,b,c, the merged metadata value should be the column name + subtype
        # e.g. for subtype 1.1, a=a1.1; b=b1.1; c=c1.1
        # if there was no subtype result, then the metadata field value should be NaN/None
        for column in list('abc'):
            if row['subtype'] != '':
                assert row[column] == column + row['subtype']
            else:
                assert np.isnan(row[column])
