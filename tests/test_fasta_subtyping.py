import pytest
from pandas import DataFrame, Series
import numpy as np
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_contigs_blastn, subtype_contigs_ac
from bio_hansel.const import SCHEME_FASTAS

from . import check_subtype_attrs, check_df_fasta_cols

genome_name = 'test'
scheme = 'heidelberg'


@pytest.fixture()
def fasta():
    return 'tests/data/SRR1002850_SMALL.fasta'


@pytest.fixture()
def fasta_gz():
    return 'tests/data/SRR1002850_SMALL.fasta.gz'


@pytest.fixture()
def exp_subtype(fasta):
    return Subtype(scheme=scheme,
                   scheme_version=SCHEME_FASTAS[scheme]['version'],
                   sample=genome_name,
                   file_path=fasta,
                   subtype='2.2.2.2.1.4',
                   inconsistent_subtypes=None,
                   n_tiles_matching_all=202,
                   n_tiles_matching_all_expected='202',
                   n_tiles_matching_positive=17,
                   n_tiles_matching_positive_expected='17',
                   n_tiles_matching_subtype=3,
                   n_tiles_matching_subtype_expected='3',
                   qc_status='PASS')


def test_heidelberg_fasta_ac(fasta, fasta_gz, exp_subtype):
    st, df = subtype_contigs_ac(scheme=scheme, input_fasta=fasta, genome_name=genome_name)
    stgz, dfgz = subtype_contigs_ac(scheme=scheme, input_fasta=fasta_gz, genome_name=genome_name)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_heidelberg_fasta_blastn(fasta, fasta_gz, exp_subtype):
    st, df = subtype_contigs_blastn(fasta_path=fasta, genome_name=genome_name, scheme=scheme)
    stgz, dfgz = subtype_contigs_blastn(fasta_path=fasta_gz, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, exp_subtype)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)
