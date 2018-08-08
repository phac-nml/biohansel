import pandas as pd
from pandas.util.testing import assert_frame_equal
from biohansel.create.read_vcf import read_vcf

"""Get the vcf file and output that proper pandas dataframe with the 2-state SNVs extracted
"""
test_vcfpath = "tests/data/create/test.vcf"
expected_sequence_df = pd.read_csv("tests/data/create/expected_sequence_df.csv", index_col=0)
expected_binary_df = pd.read_csv("tests/data/create/expected_binary_df.csv", index_col=0)


def test_readvcf():
    test_sequence_df, test_binary_df = read_vcf(test_vcfpath)
    assert_frame_equal(expected_binary_df, test_binary_df)
    assert_frame_equal(expected_sequence_df, test_sequence_df)
