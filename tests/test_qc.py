# -*- coding: utf-8 -*-

from pandas import DataFrame

from bio_hansel.qc.const import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads_ac

genome_name = 'test'


def test_intermediate_subtype():
    scheme = 'enteritidis'
    fastq = 'tests/data/Retro1000data/10-1358.fastq'
    st, df = subtype_reads_ac(reads=fastq, genome_name=genome_name, scheme='enteritidis')
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert QC.INTERMEDIATE_SUBTYPE_WARNING in st.qc_message
    assert "Total subtype matches observed (n=3) vs expected (n=6)" in st.qc_message
    assert st.qc_status == QC.WARNING


def test_missing_tiles():
    scheme = 'heidelberg'
    fastq = 'tests/data/SRR1696752/SRR1696752.fastq'
    st, df = subtype_reads_ac(reads=fastq, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.is_fastq_input()
    assert st.scheme == scheme
    assert QC.MISSING_TILES_ERROR_1 in st.qc_message
    assert 'Low coverage depth (10.9 < 20.0 expected)' in st.qc_message
    assert st.qc_status == QC.FAIL


def test_mixed_tiles():
    scheme = 'heidelberg'
    fastqs = ['tests/data/SRR3392166/SRR3392166.fastq', 'tests/data/SRR3392166/SRR3392166.fastq']
    st, df = subtype_reads_ac(reads=fastqs, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert QC.MIXED_SAMPLE_ERROR_2 in st.qc_message
    assert 'Mixed subtypes found: "1; 2; 2.1"' in st.qc_message
    assert st.qc_status == QC.FAIL
