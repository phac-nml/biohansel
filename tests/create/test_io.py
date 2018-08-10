
import pandas as pd

from Bio import SeqIO
from pandas.util.testing import assert_frame_equal

from biohansel.create.io.parsers import parse_vcf, parse_sequence_file

"""Get the vcf file and output that proper pandas dataframe with the 2-state SNVs extracted
"""
test_vcfpath = "tests/data/create/test.vcf"
expected_sequence_df = pd.read_csv("tests/data/create/expected_sequence_df.csv", index_col=0)
expected_binary_df = pd.read_csv("tests/data/create/expected_binary_df.csv", index_col=0)


def test_parse_vcf():
    test_sequence_df, test_binary_df = parse_vcf(test_vcfpath)
    assert_frame_equal(expected_binary_df, test_binary_df)
    assert_frame_equal(expected_sequence_df, test_sequence_df)

def test_parse_sequence_file():
    test_records = [
        record
        for record in SeqIO.parse("tests/data/create/sequence_ecoli.gb", "genbank")
    ]
    test_record_dict = {}
    test_record_dict[test_records[0].name] = test_records[0].seq

    test_result = parse_sequence_file("tests/data/create/sequence_ecoli.gb", "genbank")
    assert (test_result == test_record_dict)


