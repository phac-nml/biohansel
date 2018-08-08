import pandas as pd

from Bio import SeqIO

from biohansel.create.write_sequence import get_sequences, get_subsequences, read_sequence_file, \
    get_sequence_string


def test_get_sequences():
    test_sequence_df = pd.read_csv("tests/data/create/sequence_df.csv", index_col=0)
    test_length = 16
    test_list = []
    test_chrom_df = pd.read_csv("tests/data/create/curr_chrom_df.csv", index_col=0)

    test_record_dict = {}
    for record in SeqIO.parse("tests/data/create/sequence_ecoli.gb", "genbank"):
        test_record_dict[record.name] = record.seq

    test_result = get_sequences(test_sequence_df, test_length, test_record_dict)
    test_list.append(test_chrom_df)

    for test_items, expected_items in zip(test_list, test_result):
        pd.util.testing.assert_frame_equal(test_items, expected_items)


# def test_write_sequences():


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


def test_read_sequence_file():
    test_records = [
        record
        for record in SeqIO.parse("tests/data/create/sequence_ecoli.gb", "genbank")
    ]
    test_record_dict = {}
    test_record_dict[test_records[0].name] = test_records[0].seq

    test_result = read_sequence_file("tests/data/create/sequence_ecoli.gb", "genbank")
    assert (test_result == test_record_dict)


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
