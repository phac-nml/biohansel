import pytest
from pandas import DataFrame

from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads_ac
from bio_hansel.const import SCHEME_FASTAS
from bio_hansel.qc.const import QC

from . import check_df_fastq_cols, check_subtype_attrs

genome_name = 'test'
scheme_heidelberg = 'heidelberg'
scheme_enteritidis = 'enteritidis'

fastq_heidelberg_pass = 'tests/data/SRR5646583_SMALL.fastq'
fastq_gz_heidelberg_pass = 'tests/data/SRR5646583_SMALL.fastq.gz'

fastqs_enteritidis_fail = ['tests/data/inconsistent_reads_fwd.fastq', 'tests/data/inconsistent_reads_rvs.fastq']


@pytest.fixture()
def subtype_heidelberg_pass():
    return Subtype(scheme=scheme_heidelberg,
                   scheme_version=SCHEME_FASTAS[scheme_heidelberg]['version'],
                   sample=genome_name,
                   subtype='2.2.1.1.1.1',
                   file_path=fastq_heidelberg_pass,
                   are_subtypes_consistent=True,
                   n_tiles_matching_all=202,
                   n_tiles_matching_all_expected='202',
                   n_tiles_matching_positive=20,
                   n_tiles_matching_positive_expected='20',
                   n_tiles_matching_subtype=2,
                   n_tiles_matching_subtype_expected='2',
                   qc_status=QC.PASS)


@pytest.fixture()
def subtype_enteritidis_fail():
    return Subtype(scheme=scheme_enteritidis,
                   scheme_version=SCHEME_FASTAS[scheme_enteritidis]['version'],
                   sample=genome_name,
                   subtype='2.2.4.1',
                   file_path=fastq_heidelberg_pass,
                   are_subtypes_consistent=False,
                   n_tiles_matching_all=182,
                   n_tiles_matching_all_expected='224',
                   n_tiles_matching_positive=22,
                   n_tiles_matching_positive_expected='21',
                   n_tiles_matching_subtype=5,
                   n_tiles_matching_subtype_expected='6',
                   qc_status=QC.FAIL)


def test_heidelberg_scheme_vs_qc_passing_reads_with_ac(subtype_heidelberg_pass):
    st, df = subtype_reads_ac(reads=fastq_heidelberg_pass, genome_name=genome_name, scheme=scheme_heidelberg)
    stgz, dfgz = subtype_reads_ac(reads=fastq_gz_heidelberg_pass, genome_name=genome_name, scheme=scheme_heidelberg)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_heidelberg_pass)
    check_df_fastq_cols(df)
    check_df_fastq_cols(dfgz)
