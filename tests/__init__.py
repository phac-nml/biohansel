import numpy as np

attrs_to_check = """
scheme
scheme_version
sample
subtype
are_subtypes_consistent
n_kmers_matching_all
n_kmers_matching_all_expected
n_kmers_matching_positive
n_kmers_matching_positive_expected
n_kmers_matching_subtype
n_kmers_matching_subtype_expected
qc_status
""".strip().split('\n')

exp_fasta_cols = ['kmername', 'contig_id', 'refposition', 'seq', 'subtype', 'match_index', 'is_revcomp',
                  'is_pos_kmer', 'sample', 'file_path', 'scheme', 'scheme_version', 'qc_status', 'qc_message']

exp_fastq_cols = ['kmername', 'refposition', 'subtype', 'seq', 'freq', 'is_pos_kmer', 'is_kmer_freq_okay',
                  'kmer_fraction','total_refposition_kmer_frequency','is_kmer_fraction_okay',
                  'sample', 'file_path', 'scheme', 'scheme_version', 'qc_status', 'qc_message']


def check_subtype_attrs(*sts):
    for a in attrs_to_check:
        assert len({st.__getattribute__(a) for st in sts}) == 1, \
            'Mismatching values for attr "{}": "{}"\n{}'.format(a,
                                                                [st.__getattribute__(a) for st in sts],
                                                                sts)


def check_df_fasta_cols(df):
    assert np.all(df.columns.isin(exp_fasta_cols))


def check_df_fastq_cols(df):
    assert np.all(df.columns.isin(exp_fastq_cols))
