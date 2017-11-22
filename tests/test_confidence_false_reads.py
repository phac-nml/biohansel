import pytest
from pandas import DataFrame, Series
import numpy as np

from bio_hansel.quality_check.const import FAIL_MESSAGE
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_reads
from bio_hansel.subtyping_params import SubtypingParams
from bio_hansel.utils import SCHEME_FASTAS


@pytest.fixture()
def test_genomes():
    return ['tests/data/inconsistent_reads_fwd.fastq', 'tests/data/inconsistent_reads_rvs.fastq']


def test_confidence_false(test_genomes):
    genome_name = 'test'
    scheme = 'enteritidis'
    st, df = subtype_reads(scheme='enteritidis', reads=test_genomes, genome_name=genome_name, subtyping_params=SubtypingParams())
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)

    assert st.scheme == scheme
    assert st.scheme_version == SCHEME_FASTAS[scheme]['version']
    assert st.sample == genome_name
    assert st.subtype == '2.2.2.4.1'
    assert st.are_subtypes_consistent is False
    assert st.n_tiles_matching_all == 208
    assert st.n_tiles_matching_all_expected == '188'
    assert st.n_tiles_matching_positive == 26
    assert st.n_tiles_matching_positive_expected == '25'
    assert st.n_tiles_matching_subtype == 5
    assert st.n_tiles_matching_subtype_expected == '6'
    assert st.qc_status == FAIL_MESSAGE

    exp_cols = ['tilename', 'freq', 'refposition', 'subtype',
                'is_pos_tile', 'is_kmer_freq_okay', 'sample', 'file_path', 'scheme', 'scheme_version',
                'qc_status', 'qc_message']
    df_cols = df.columns  # type: Series
    assert np.all(df_cols.isin(exp_cols))
