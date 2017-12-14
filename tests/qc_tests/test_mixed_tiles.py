import pytest
from pandas import DataFrame

from bio_hansel.quality_check.const import MIXED_SAMPLE_ERROR_2
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads_jellyfish
from bio_hansel.subtyping_params import SubtypingParams


@pytest.fixture()
def test_genomes():
    return ['tests/data/SRR3392166/SRR3392166.fastq', 'tests/data/SRR3392166/SRR3392166.fastq']


def test_mixed_tiles(test_genomes):
    genome_name = 'test'
    scheme = 'heidelberg'
    st, df = subtype_reads_jellyfish(reads=test_genomes, genome_name=genome_name, scheme='heidelberg',
                                     subtyping_params=SubtypingParams())
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert MIXED_SAMPLE_ERROR_2 in st.qc_message
    assert "Mixed subtypes found: ['1', '2', '2.1']" in st.qc_message
    assert st.qc_status == 'FAIL'
