import pytest
from pandas import DataFrame

from bio_hansel.const import SCHEME_FASTAS
from bio_hansel.qc.const import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads,subtype_contigs
from . import check_df_fastq_cols, check_subtype_attrs

genome_name = 'test'
scheme_heidelberg = 'heidelberg'
scheme_enteritidis = 'enteritidis'
scheme_typhi = 'typhi'
scheme_tuberculosis = 'tb_lineage'
scheme_typhimurium = 'typhimurium'

fasta_heidelberg_pass = 'tests/data/SRR1002850_lowercase.fasta'
fastq_heidelberg_pass = 'tests/data/SRR5646583_SMALL.fastq'
fastq_gz_heidelberg_pass = 'tests/data/SRR5646583_SMALL.fastq.gz'

fastqs_enteritidis_fail = ['tests/data/inconsistent_reads_fwd.fastq', 'tests/data/inconsistent_reads_rvs.fastq']

fasta_typhi_pass = 'tests/data/AE014613.1.fasta'

fasta_tb_pass = 'tests/data/AP018036.1.fasta'

fasta_typhimurium_pass = 'tests/data/typhimurium2.2.3.3.fasta'


@pytest.fixture()
def subtype_heidelberg_pass():
    return Subtype(scheme=scheme_heidelberg,
                   scheme_version=SCHEME_FASTAS[scheme_heidelberg]['version'],
                   sample=genome_name,
                   subtype='2.2.1.1.1.1',
                   file_path=fastq_heidelberg_pass,
                   are_subtypes_consistent=True,
                   n_kmers_matching_all=202,
                   n_kmers_matching_all_expected='202',
                   n_kmers_matching_positive=20,
                   n_kmers_matching_positive_expected='20',
                   n_kmers_matching_subtype=2,
                   n_kmers_matching_subtype_expected='2',
                   qc_status=QC.PASS)


@pytest.fixture()
def subtype_enteritidis_fail():
    return Subtype(scheme=scheme_enteritidis,
                   scheme_version=SCHEME_FASTAS[scheme_enteritidis]['version'],
                   sample=genome_name,
                   subtype='2.2.4.1',
                   file_path=fastq_heidelberg_pass,
                   are_subtypes_consistent=False,
                   n_kmers_matching_all=182,
                   n_kmers_matching_all_expected='224',
                   n_kmers_matching_positive=22,
                   n_kmers_matching_positive_expected='21',
                   n_kmers_matching_subtype=5,
                   n_kmers_matching_subtype_expected='6',
                   qc_status=QC.FAIL)

@pytest.fixture()
def subtype_heidelberg_SRR1002850_pass():
    return Subtype(scheme=scheme_heidelberg,
                   scheme_version=SCHEME_FASTAS[scheme_heidelberg]['version'],
                   sample=genome_name,
                   subtype='2.2.2.2.1.4',
                   file_path=fastq_heidelberg_pass,
                   are_subtypes_consistent=True,
                   n_kmers_matching_all=202,
                   n_kmers_matching_all_expected='202',
                   n_kmers_matching_positive=17,
                   n_kmers_matching_positive_expected='17',
                   n_kmers_matching_subtype=3,
                   n_kmers_matching_subtype_expected='3',
                   qc_status=QC.PASS)

@pytest.fixture()
def subtype_typhi_AE014613_pass():
    return Subtype(scheme=scheme_typhi,
                   scheme_version=SCHEME_FASTAS[scheme_typhi]['version'],
                   sample=genome_name,
                   subtype='2.3.6.1',
                   file_path=fasta_typhi_pass,
                   are_subtypes_consistent=True,
                   n_kmers_matching_all=68,
                   n_kmers_matching_all_expected='68',
                   n_kmers_matching_positive=4,
                   n_kmers_matching_positive_expected='4',
                   n_kmers_matching_subtype=1,
                   n_kmers_matching_subtype_expected='1',
                   qc_status=QC.PASS)

@pytest.fixture()
def subtype_tb_AP018036_pass():
    return Subtype(scheme=scheme_tuberculosis,
                   scheme_version=SCHEME_FASTAS[scheme_tuberculosis]['version'],
                   sample=genome_name,
                   subtype='2.2.1',
                   file_path=fasta_tb_pass,
                   are_subtypes_consistent=True,
                   n_kmers_matching_all=62,
                   n_kmers_matching_all_expected='62',
                   n_kmers_matching_positive=3,
                   n_kmers_matching_positive_expected='3',
                   n_kmers_matching_subtype=1,
                   n_kmers_matching_subtype_expected='1',
                   qc_status=QC.PASS)

@pytest.fixture()
def subtype_typhimurium_pass():
    return Subtype(scheme=scheme_typhimurium,
                   scheme_version=SCHEME_FASTAS[scheme_typhimurium]['version'],
                   sample=genome_name,
                   subtype='2.2.3.3',
                   file_path=fasta_typhimurium_pass,
                   are_subtypes_consistent=True,
                   n_kmers_matching_all=429,
                   n_kmers_matching_all_expected='430',
                   n_kmers_matching_positive=19,
                   n_kmers_matching_positive_expected='19',
                   n_kmers_matching_subtype=5,
                   n_kmers_matching_subtype_expected='5',
                   qc_status=QC.PASS)


def test_heidelberg_scheme_vs_qc_passing_reads_with_ac(subtype_heidelberg_pass):
    st, df = subtype_reads(reads=fastq_heidelberg_pass, genome_name=genome_name, scheme=scheme_heidelberg)
    stgz, dfgz = subtype_reads(reads=fastq_gz_heidelberg_pass, genome_name=genome_name, scheme=scheme_heidelberg)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_heidelberg_pass)
    check_df_fastq_cols(df)
    check_df_fastq_cols(dfgz)



def test_heidelberg_scheme_and_lowcase_seq_inputs(subtype_heidelberg_SRR1002850_pass):
    st, df = subtype_contigs(fasta_path=fasta_heidelberg_pass, genome_name=genome_name, scheme=scheme_heidelberg)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    check_subtype_attrs(st,subtype_heidelberg_SRR1002850_pass)

def test_typhi_scheme(subtype_typhi_AE014613_pass):
    st, df = subtype_contigs(fasta_path=fasta_typhi_pass, genome_name=genome_name, scheme=scheme_typhi)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    check_subtype_attrs(st,subtype_typhi_AE014613_pass)

def test_tuberculosis_scheme(subtype_tb_AP018036_pass):
    st, df = subtype_contigs(fasta_path=fasta_tb_pass, genome_name=genome_name, scheme=scheme_tuberculosis)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    check_subtype_attrs(st,subtype_tb_AP018036_pass)

def test_typhimurium_scheme(subtype_typhimurium_pass):
    st, df = subtype_contigs(fasta_path=fasta_typhimurium_pass, genome_name=genome_name, scheme=scheme_typhimurium)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    check_subtype_attrs(st,subtype_typhimurium_pass)