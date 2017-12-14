import pytest
from pandas import DataFrame

from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads_jellyfish, subtype_reads_ac
from bio_hansel.const import SCHEME_FASTAS

from . import check_df_fastq_cols, check_subtype_attrs

scheme = 'heidelberg'
genome_name = 'test'


@pytest.fixture
def fastq():
    return 'tests/data/SRR5646583_SMALL.fastq'


@pytest.fixture()
def fastq_gz():
    return 'tests/data/SRR5646583_SMALL.fastq.gz'

@pytest.fixture()
def exp_subtype(fastq):
    return Subtype(scheme=scheme,
                   scheme_version=SCHEME_FASTAS[scheme]['version'],
                   sample=genome_name,
                   subtype='2.2.1.1.1.1',
                   file_path=fastq,
                   are_subtypes_consistent=True,
                   n_tiles_matching_all=202,
                   n_tiles_matching_all_expected='202',
                   n_tiles_matching_positive=20,
                   n_tiles_matching_positive_expected='20',
                   n_tiles_matching_subtype=2,
                   n_tiles_matching_subtype_expected='2',
                   qc_status='PASS')



def test_heidelberg_fastq_jellyfish(fastq, fastq_gz, exp_subtype):
    st, df = subtype_reads_jellyfish(reads=fastq, genome_name=genome_name, scheme=scheme, threads=4)
    stgz, dfgz = subtype_reads_jellyfish(reads=fastq_gz, genome_name=genome_name, scheme=scheme, threads=4)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fastq_cols(df)
    check_df_fastq_cols(dfgz)

def test_heidelberg_fastq_ac(fastq, fastq_gz, exp_subtype):
    st, df = subtype_reads_ac(reads=fastq, genome_name=genome_name, scheme=scheme)
    stgz, dfgz = subtype_reads_ac(reads=fastq_gz, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fastq_cols(df)
    check_df_fastq_cols(dfgz)
