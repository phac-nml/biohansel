import pytest
from pandas import DataFrame, Series
import numpy as np
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_contigs_blastn, subtype_contigs_ac
from bio_hansel.const import SCHEME_FASTAS
from bio_hansel.qc.const import QC

from . import check_subtype_attrs, check_df_fasta_cols

genome_name = 'test'
scheme_heidelberg = 'heidelberg'
scheme_enteritidis = 'enteritidis'

fasta_heidelberg_pass = 'tests/data/SRR1002850_SMALL.fasta'
fasta_gz_heidelberg_pass = 'tests/data/SRR1002850_SMALL.fasta.gz'

# input contigs that should give an unconfident result and a QC fail
fasta_enteritidis_unconfident = 'tests/data/fail-qc-unconfident-subtype.fasta'
fasta_gz_enteritidis_unconfident = 'tests/data/fail-qc-unconfident-subtype.fasta.gz'
# input contigs that should give an inconsistent result and a QC fail
fasta_enteritidis_fail = 'tests/data/SRR1958005.fasta'
fasta_gz_enteritidis_fail = 'tests/data/SRR1958005.fasta.gz'


@pytest.fixture()
def subtype_heidelberg_pass():
    return Subtype(scheme=scheme_heidelberg,
                   scheme_version=SCHEME_FASTAS[scheme_heidelberg]['version'],
                   sample=genome_name,
                   file_path=fasta_heidelberg_pass,
                   subtype='2.2.2.2.1.4',
                   are_subtypes_consistent=True,
                   inconsistent_subtypes=None,
                   n_tiles_matching_all=202,
                   n_tiles_matching_all_expected='202',
                   n_tiles_matching_positive=17,
                   n_tiles_matching_positive_expected='17',
                   n_tiles_matching_subtype=3,
                   n_tiles_matching_subtype_expected='3',
                   qc_status=QC.PASS)


@pytest.fixture()
def subtype_enteritidis_fail_unconfident():
    return Subtype(scheme=scheme_enteritidis,
                   scheme_version=SCHEME_FASTAS[scheme_enteritidis]['version'],
                   sample=genome_name,
                   subtype='2.1.1',
                   file_path=fasta_enteritidis_unconfident,
                   are_subtypes_consistent=True,
                   n_tiles_matching_all=176,
                   n_tiles_matching_all_expected='188',
                   n_tiles_matching_positive=9,
                   n_tiles_matching_positive_expected='9',
                   n_tiles_matching_subtype=1,
                   n_tiles_matching_subtype_expected='1',
                   qc_status=QC.FAIL)


@pytest.fixture()
def subtype_enteritidis_fail():
    return Subtype(scheme=scheme_enteritidis,
                   scheme_version=SCHEME_FASTAS[scheme_enteritidis]['version'],
                   sample=genome_name,
                   subtype='2.1; 2.2',
                   file_path=fasta_enteritidis_fail,
                   are_subtypes_consistent=False,
                   n_tiles_matching_all=188,
                   n_tiles_matching_all_expected='188;188',
                   n_tiles_matching_positive=9,
                   n_tiles_matching_positive_expected='8;10',
                   n_tiles_matching_subtype=3,
                   n_tiles_matching_subtype_expected='2;4',
                   qc_status=QC.FAIL)


def test_heidelberg_fasta_ac(subtype_heidelberg_pass):
    st, df = subtype_contigs_ac(fasta_path=fasta_heidelberg_pass,
                                genome_name=genome_name,
                                scheme=scheme_heidelberg)
    stgz, dfgz = subtype_contigs_ac(fasta_path=fasta_gz_heidelberg_pass,
                                    genome_name=genome_name,
                                    scheme=scheme_heidelberg)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_heidelberg_pass)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_heidelberg_fasta_blastn(subtype_heidelberg_pass):
    st, df = subtype_contigs_blastn(fasta_path=fasta_heidelberg_pass,
                                    genome_name=genome_name,
                                    scheme=scheme_heidelberg)
    stgz, dfgz = subtype_contigs_blastn(fasta_path=fasta_gz_heidelberg_pass,
                                        genome_name=genome_name,
                                        scheme=scheme_heidelberg)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_heidelberg_pass)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_enteritidis_scheme_vs_qc_failing_contigs_unconfident_blastn(subtype_enteritidis_fail_unconfident):
    st, df = subtype_contigs_blastn(fasta_path=fasta_enteritidis_unconfident,
                                    genome_name=genome_name,
                                    scheme=scheme_enteritidis)
    stgz, dfgz = subtype_contigs_blastn(fasta_path=fasta_gz_enteritidis_unconfident,
                                        genome_name=genome_name,
                                        scheme=scheme_enteritidis)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_enteritidis_fail_unconfident)
    assert 'Unconfident Results Error 4' in st.qc_message
    assert 'Unconfident Results Error 4' in stgz.qc_message
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_enteritidis_scheme_vs_qc_failing_contigs_unconfident_ac(subtype_enteritidis_fail_unconfident):
    st, df = subtype_contigs_ac(fasta_path=fasta_enteritidis_unconfident,
                                genome_name=genome_name,
                                scheme=scheme_enteritidis)
    stgz, dfgz = subtype_contigs_ac(fasta_path=fasta_gz_enteritidis_unconfident,
                                    genome_name=genome_name,
                                    scheme=scheme_enteritidis)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_enteritidis_fail_unconfident)
    assert 'Unconfident Results Error 4' in st.qc_message
    assert 'Unconfident Results Error 4' in stgz.qc_message
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_blastn_vs_bad_contigs(subtype_enteritidis_fail):
    st, df = subtype_contigs_blastn(fasta_path=fasta_enteritidis_fail,
                                    genome_name=genome_name,
                                    scheme=scheme_enteritidis)
    stgz, dfgz = subtype_contigs_blastn(fasta_path=fasta_gz_enteritidis_fail,
                                        genome_name=genome_name,
                                        scheme=scheme_enteritidis)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_enteritidis_fail)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)


def test_ac_vs_bad_contigs(subtype_enteritidis_fail):
    st, df = subtype_contigs_ac(fasta_path=fasta_enteritidis_fail,
                                genome_name=genome_name,
                                scheme=scheme_enteritidis)
    stgz, dfgz = subtype_contigs_ac(fasta_path=fasta_gz_enteritidis_fail,
                                    scheme=scheme_enteritidis,
                                    genome_name=genome_name)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert isinstance(stgz, Subtype)
    assert isinstance(dfgz, DataFrame)
    check_subtype_attrs(st, stgz, subtype_enteritidis_fail)
    check_df_fasta_cols(df)
    check_df_fasta_cols(dfgz)
