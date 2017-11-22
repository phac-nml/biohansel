import pytest
from pandas import DataFrame

from bio_hansel.quality_check.const import MISSING_TILES_ERROR_1A
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads
from bio_hansel.subtyping_params import SubtypingParams


@pytest.fixture()
def test_genomes():
    return ['tests/data/SRR1696752/SRR1696752_forward.fastqsanger', 'tests/data/SRR1696752/SRR1696752_reverse.fastqsanger']


def test_missing_tiles(test_genomes):
    genome_name = 'test'
    scheme = 'heidelberg'
    st, df = subtype_reads(scheme='heidelberg', reads=test_genomes, genome_name=genome_name, subtyping_params=SubtypingParams())
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert MISSING_TILES_ERROR_1A in st.qc_message
    assert 'Average calculated tile coverage = 12.5' in st.qc_message
    assert st.qc_status == 'FAIL'
