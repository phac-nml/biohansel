import pytest
from pandas import DataFrame, Series
import numpy as np

from bio_hansel.const import MIXED_SUBTYPE_ERROR, OK_NUM_TILES
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads
from bio_hansel.utils import SCHEME_FASTAS


@pytest.fixture()
def test_genomes():
    # The data to provided for false confidence is not here yet, substitute this file name with the correct data.
    return ['tests/data/fwd_inconsistent.fastqsanger', 'tests/data/fwd_inconsistent.fastqsanger']


def test_confidence_false(test_genomes):
    genome_name = 'test'
    scheme = 'enteritidis'
    st, df = subtype_reads(scheme='enteritidis', reads=test_genomes, genome_name=genome_name)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)

    assert st.scheme == scheme
    assert st.scheme_version == SCHEME_FASTAS[scheme]['version']
    assert st.sample == genome_name
    assert st.subtype == '2.1; 2.2'
    assert st.are_subtypes_consistent is False
    assert st.n_tiles_matching_all == 186
    assert st.n_tiles_matching_all_expected == '188;188'
    assert st.n_tiles_matching_positive == 9
    assert st.n_tiles_matching_positive_expected == '8;10'
    assert st.n_tiles_matching_subtype == 3
    assert st.n_tiles_matching_subtype_expected == '2;4'
    assert st.confident_is_subtype is MIXED_SUBTYPE_ERROR
    assert st.reached_min_tiles is OK_NUM_TILES
    exp_cols = ['tilename', 'freq', 'refposition', 'subtype',
                'is_pos_tile', 'is_kmer_freq_okay', 'sample', 'file_path', 'scheme', 'scheme_version',
                'reached_min_tiles', 'is_confident']
    df_cols = df.columns  # type: Series
    assert np.all(df_cols.isin(exp_cols))
