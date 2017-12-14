import numpy as np

attrs_to_check = """
scheme
scheme_version
sample
subtype
are_subtypes_consistent
n_tiles_matching_all
n_tiles_matching_all_expected
n_tiles_matching_positive
n_tiles_matching_positive_expected
n_tiles_matching_subtype
n_tiles_matching_subtype_expected
qc_status
""".strip().split('\n')


exp_fasta_cols = ['tilename', 'contig_id', 'refposition', 'subtype', 'match_index', 'is_revcomp',
            'is_pos_tile', 'sample', 'file_path', 'scheme', 'scheme_version', 'qc_status', 'qc_message']

exp_fastq_cols = ['tilename', 'contig_id', 'refposition', 'subtype', 'freq', 'is_pos_tile', 'is_kmer_freq_okay',
                  'sample', 'file_path', 'scheme', 'scheme_version', 'qc_status', 'qc_message']

def check_subtype_attrs(*sts):
    for a in attrs_to_check:
        assert len({st.__getattribute__(a) for st in sts}) == 1


def check_df_fasta_cols(df):
    assert np.all(df.columns.isin(exp_fasta_cols))


def check_df_fastq_cols(df):
    assert np.all(df.columns.isin(exp_fastq_cols))