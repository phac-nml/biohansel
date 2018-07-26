import pandas as pd
from pandas.util.testing import assert_frame_equal
from bio_hansel.scripts.readvcf import read_vcf
"""Get the vcf file and output that proper pandas dataframe with the 2-state SNVs extracted
"""
test_vcfpath = "tests/data/test.vcf"
expected_df = pd.read_csv("tests/data/expected.csv")


def test_readvcf():

    test_df = read_vcf(test_vcfpath)
    pd.testing.assert_frame_equal(expected_df, test_df)
