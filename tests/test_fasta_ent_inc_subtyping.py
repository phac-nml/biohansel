import numpy as np
import pytest
from pandas import DataFrame, Series

from bio_hansel.const import SCHEME_FASTAS
from bio_hansel.quality_check import FAIL_MESSAGE
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_contigs_blastn, subtype_contigs_ac
from . import check_subtype_attrs, check_df_fasta_cols

scheme = 'enteritidis'
genome_name = 'test'


@pytest.fixture()
def fasta():
    return 'tests/data/SRR1958005.fasta'


@pytest.fixture()
def fasta_gz():
    return 'tests/data/SRR1958005.fasta.gz'


@pytest.fixture()
def exp_subtype(fasta):
    return Subtype(scheme=scheme,
                   scheme_version=SCHEME_FASTAS[scheme]['version'],
                   sample=genome_name,
                   subtype='2.1; 2.2',
                   file_path=fasta,
                   are_subtypes_consistent=False,
                   n_tiles_matching_all=188,
                   n_tiles_matching_all_expected='188;188',
                   n_tiles_matching_positive=9,
                   n_tiles_matching_positive_expected='8;10',
                   n_tiles_matching_subtype=3,
                   n_tiles_matching_subtype_expected='2;4',
                   qc_status=FAIL_MESSAGE)


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
    st, df = subtype_contigs_ac(scheme=scheme, input_fasta=fasta, genome_name=genome_name)
    stgz, dfgz = subtype_contigs_ac(input_fasta=fasta_gz, scheme=scheme, genome_name=genome_name)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)
