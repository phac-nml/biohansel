import pytest
import pandas as pd
from biohansel.create.group_snvs import group_snvs
import numpy.testing as npt

def test_group_snvs():
    """
    Tests whether or not the group_snvs will provide the expected groupings
    """
    test_dict = {'mysnpsSRR6683541': '1', 'mysnpsSRR6683736': '1', 'mysnpsSRR6683914': '1', 'Reference': '1', 'mysnpsSRR6683916': '2'}
    test_binary_df = pd.read_csv('tests/data/create/expected_binary_df.csv', index_col=0)
    test_sequence_df = pd.read_csv('tests/data/create/expected_sequence_df.csv', index_col=0)
    group1=pd.read_csv("tests/data/create/group1.csv", dtype={'CHROM':str, 'REF':str, 'ALT':str, 'mysnpsSRR6683541':float} ,index_col=0)
    group2=pd.read_csv("tests/data/create/group2.csv", dtype={'CHROM':str, 'REF':str, 'ALT':str, 'mysnpsSRR6683916':float}, index_col=0)

    
    test_result=group_snvs(test_binary_df, test_sequence_df, test_dict)

    pd.testing.assert_frame_equal(group1, test_result["1"])
    pd.testing.assert_frame_equal(group2, test_result["2"]) 