import pytest
from pandas import DataFrame
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_contigs_blastn, subtype_contigs_ac
from bio_hansel.const import SCHEME_FASTAS

from . import check_subtype_attrs, check_df_fasta_cols

genome_name = 'test'
scheme = 'enteritidis'


@pytest.fixture()
def fasta():
    return 'tests/data/SRR6126859.fasta'


@pytest.fixture()
def fasta_gz():
    return 'tests/data/SRR6126859.fasta.gz'


@pytest.fixture()
def exp_subtype(fasta):
    return Subtype(scheme=scheme,
                   scheme_version=SCHEME_FASTAS[scheme]['version'],
                   sample=genome_name,
                   subtype='2.2.2',
                   file_path=fasta,
                   are_subtypes_consistent=True,
                   n_tiles_matching_all=188,
                   n_tiles_matching_all_expected='188',
                   n_tiles_matching_positive=16,
                   n_tiles_matching_positive_expected='16',
                   n_tiles_matching_subtype=6,
                   n_tiles_matching_subtype_expected='6',
                   qc_status='PASS')


def test_ent_fasta_subtyping_blastn(fasta, fasta_gz, exp_subtype):
    st, df = subtype_contigs_blastn(fasta_path=fasta, genome_name=genome_name, scheme=scheme)
    stgz, dfgz = subtype_contigs_blastn(fasta_path=fasta_gz, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_ent_fasta_subtyping_ac(fasta, fasta_gz, exp_subtype):
    st, df = subtype_contigs_ac(input_fasta=fasta, genome_name=genome_name, scheme=scheme)
    stgz, dfgz = subtype_contigs_ac(input_fasta=fasta_gz, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)
