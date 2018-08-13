import pandas as pd

from Bio import SeqIO

from biohansel.create.schema_generator import get_sequences, get_subsequences,get_sequence_string, group_snvs


def test_get_sequences():
    test_sequence_df = pd.read_csv("tests/data/create/sequence_df.csv", index_col=0)
    test_length = 32
    test_list = []
    test_chrom_df = pd.read_csv("tests/data/create/curr_chrom_df.csv", index_col=0)

    test_record_dict = {}
    for record in SeqIO.parse("tests/data/create/sequence_ecoli.gb", "genbank"):
        test_record_dict[record.name] = record.seq

    test_result = get_sequences(test_sequence_df, test_length, test_record_dict)
    test_list.append(test_chrom_df)

    for test_items, expected_items in zip(test_list, test_result):
        pd.util.testing.assert_frame_equal(test_items, expected_items)

def test_get_subsequences():
    test_records = [
        record
        for record in SeqIO.parse("tests/data/create/sequence_ecoli.gb", "genbank")
    ]
    expected_result = "CAGAACGTTTTCTGCGGGTTGCCGATATTCTGG"
    test_seq = test_records[0].seq
    test_length = 16
    test_position = 410
    test_max_sequence_value = len(test_seq)
    test_result = get_subsequences(test_position, test_seq, test_length, test_max_sequence_value)

    assert (test_result == expected_result)

def test_get_sequence_string():
    test_string = """>(NC_002695)576-2
TATTTTTGCCGAACTTTTGACGGGACTCGCCGC
>negative(NC_002695)576-2
TATTTTTGCCGAACTTCTGACGGGACTCGCCGC\n"""
    test_ratio_value = 0
    test_chromosome = "NC_002695"
    test_position = 576
    test_group = 2
    test_reference_snv = "TATTTTTGCCGAACTTTTGACGGGACTCGCCGC"
    test_alternate_snv = "TATTTTTGCCGAACTTCTGACGGGACTCGCCGC"
    test_result = get_sequence_string(test_ratio_value, test_chromosome, test_position, test_group, test_reference_snv,
                                      test_alternate_snv)

    assert (test_string == test_result)

def test_group_snvs():
    """
    Tests whether or not the group_snvs will provide the expected groupings
    """
    test_dict = {'SRR6683541': '1', 'SRR6683736': '1', 'SRR6683914': '1', 'Reference': '1', 'SRR6683916': '2'}
    test_binary_df = pd.read_csv('tests/data/create/expected_binary_df.csv', index_col=0)
    test_sequence_df = pd.read_csv('tests/data/create/expected_sequence_df.csv', index_col=0)
    group1 = pd.read_csv("tests/data/create/group1.csv",
                         dtype={'CHROM': str, 'REF': str, 'ALT': str, 'SRR6683541': float}, index_col=0)
    group2 = pd.read_csv("tests/data/create/group2.csv",
                         dtype={'CHROM': str, 'REF': str, 'ALT': str, 'SRR6683916': float}, index_col=0)

    test_result = group_snvs(test_binary_df, test_sequence_df, test_dict)

    pd.testing.assert_frame_equal(group1, test_result["1"])
    pd.testing.assert_frame_equal(group2, test_result["2"])