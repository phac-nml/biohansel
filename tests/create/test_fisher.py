import pytest
import pandas as pd
from bio_hansel.scripts.fisher_test import fisherTest
import numpy.testing as npt

test_dict = {
    "SRR123455": "1",
    "SRR12345": "2",
    "SRR123212": "3",
    "SRR32323": "1",
    "SRR6969": "2"
}
test_data = pd.read_csv('tests/data/sampleData.csv')
test_data1 = pd.read_csv('tests/data/sampleData.csv')
group_1 = pd.read_csv('tests/data/group1.csv')
group_2 = pd.read_csv('tests/data/group2.csv')
group_3 = pd.read_csv('tests/data/group3.csv')


def test_fishertest():
    """
    Tests whether or not the fishertest will provide the expected groupings
    """
    result = fisherTest(test_data, test_dict)

    npt.assert_array_equal(group_1.values, result["1"].values)
    npt.assert_array_equal(group_2.values, result["2"].values)
    npt.assert_array_equal(group_3.values, result["3"].values)