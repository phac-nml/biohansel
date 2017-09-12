import pytest
from pandas import DataFrame, Series
import numpy as np

from bio_hansel.subtype import Subtype
from bio_hansel.subtyper import subtype_fasta


@pytest.fixture
def test_genome():
    return 'tests/data/FSIS1700719.fasta'


def test_fasta_subtyping(test_genome):
    genome_name = 'test'
    scheme = 'heidelberg'
    st, df = subtype_fasta(scheme='heidelberg', fasta_path=test_genome, genome_name=genome_name)
    assert isinstance(st, Subtype)
    assert isinstance(df, DataFrame)

    assert st.scheme == scheme
    assert st.sample == genome_name
    assert st.subtype == '2.2.1.1.1.1'
    assert st.inconsistent_subtypes is None
    assert st.n_tiles_matching_all == 202
    assert st.n_tiles_matching_all_total == '202'

    exp_cols = ['tilename', 'stitle', 'pident', 'length', 'mismatch', 'gapopen',
       'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen',
       'slen', 'seq', 'coverage', 'is_trunc', 'refposition', 'subtype',
       'is_pos_tile', 'sample', 'file_path', 'scheme']
    df_cols = df.columns # type: Series
    assert np.all(df_cols.isin(exp_cols))
