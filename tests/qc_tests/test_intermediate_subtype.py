import pytest
from pandas import DataFrame

from bio_hansel.quality_check.const import INTERMEDIATE_SUBTYPE_WARNING, AMBIGUOUS_RESULTS_ERROR_3
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads
from bio_hansel.subtyping_params import SubtypingParams


@pytest.fixture()
def test_genomes():
    return ['tests/data/Retro1000data/10-1358.fastq', 'tests/data/Retro1000data/10-1358.fastq']


def test_intermediate_subtype(test_genomes):
    genome_name = 'test'
    scheme = 'enteritidis'
    st, df = subtype_reads(scheme='enteritidis', reads=test_genomes, genome_name=genome_name, subtyping_params=SubtypingParams())
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert INTERMEDIATE_SUBTYPE_WARNING in st.qc_message
    assert "Total subtype hits: 3 | Total subtype expected: 6" in st.qc_message
    assert st.qc_status == 'WARNING'
