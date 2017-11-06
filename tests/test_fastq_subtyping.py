import pytest
from pandas import DataFrame
import numpy as np

from bio_hansel.const import PASS_MESSAGE
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads
from bio_hansel.utils import SCHEME_FASTAS


@pytest.fixture
def test_genome():
    return 'tests/data/SRR5646583_SMALL.fastq'


def test_fastq_subtyping(test_genome):
    genome_name = 'test'
    scheme = 'heidelberg'
    st, df = subtype_reads(scheme='heidelberg', reads=test_genome, genome_name=genome_name, threads=4)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)

    assert st.scheme == scheme
    assert st.scheme_version == SCHEME_FASTAS[scheme]['version']
    assert st.sample == genome_name
    assert st.subtype == '2.2.1.1.1.1'
    assert st.are_subtypes_consistent is True
    assert st.inconsistent_subtypes is None
    assert st.n_tiles_matching_all == 202
    assert st.n_tiles_matching_all_expected == '202'
    assert st.n_tiles_matching_positive == 20
    assert st.n_tiles_matching_positive_expected == '20'
    assert st.n_tiles_matching_subtype == 2
    assert st.n_tiles_matching_subtype_expected == '2'
    assert st.qc_status == PASS_MESSAGE
    assert len(st.qc_message) == 0

    exp_cols = ['tilename', 'freq', 'refposition', 'subtype',
                'is_pos_tile', 'is_kmer_freq_okay', 'sample', 'file_path', 'scheme', 'scheme_version',
                'qc_status', 'qc_message']
    df_cols = df.columns  # type: Series
    assert np.all(df_cols.isin(exp_cols))
