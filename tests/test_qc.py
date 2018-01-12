# -*- coding: utf-8 -*-

from pandas import DataFrame

from bio_hansel.qc.const import QC
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads_ac, subtype_contigs_ac

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


def test_mixed_subtype_positive_negative_tiles_same_target():
    scheme = 'heidelberg'
    fasta = 'tests/data/fail-qc-mixed-subtype-pos-neg-tiles.fasta'
    st, df = subtype_contigs_ac(fasta_path=fasta, genome_name=genome_name, scheme=scheme)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert st.qc_status == QC.FAIL
    assert QC.MIXED_SAMPLE_ERROR_2 in st.qc_message
    expected_qc_msg = 'FAIL: Mixed Sample Error 2: Mixed subtype detected. ' \
                      'Positive and negative tiles detected for the same ' \
                      'target site ' \
                      '"202001; 600783; 1049933; 1193219; 2778621; 2904061; ' \
                      '3278067; 3867228; 4499501; 4579224; 4738855; 202001; ' \
                      '600783; 1049933; 1193219; 2778621; 2904061; 3278067; ' \
                      '3867228; 4499501; 4579224; 4738855" for subtype "1.1".'
    assert expected_qc_msg in st.qc_message
