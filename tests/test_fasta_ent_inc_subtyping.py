import pytest
from pandas import DataFrame, Series
import numpy as np

from bio_hansel.quality_check import FAIL_MESSAGE
from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_fasta
from bio_hansel.subtyping_params import SubtypingParams
from bio_hansel.const import SCHEME_FASTAS


@pytest.fixture()
def test_genome():
    return 'tests/data/SRR1958005.fasta'


def test_ent_fasta_subtyping(test_genome):
    genome_name = 'test'
    scheme = 'enteritidis'
    st, df = subtype_fasta(scheme='enteritidis', fasta_path=test_genome, genome_name=genome_name, subtyping_params=SubtypingParams())
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)
    assert st.scheme == scheme
    assert st.scheme_version == SCHEME_FASTAS[scheme]['version']
    assert st.sample == genome_name
    assert st.subtype == '2.1; 2.2'
    assert st.are_subtypes_consistent == False
    assert st.n_tiles_matching_all == 188
    assert st.n_tiles_matching_all_expected == '188;188'
    assert st.n_tiles_matching_positive == 9
    assert st.n_tiles_matching_positive_expected == '8;10'
    assert st.n_tiles_matching_subtype == 3
    assert st.n_tiles_matching_subtype_expected == '2;4'
    assert st.qc_status == FAIL_MESSAGE

    exp_cols = ['tilename', 'stitle', 'refposition', 'subtype',
       'is_pos_tile', 'sample', 'file_path', 'scheme', 'scheme_version', 'qc_status', 'qc_message']
    df_cols = df.columns # type: Series
    assert np.all(df_cols.isin(exp_cols))
