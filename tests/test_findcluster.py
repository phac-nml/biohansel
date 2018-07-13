import pytest
import pandas as pd
from bio_hansel.scripts.findcluster import find_clusters
from bio_hansel.scripts.readvcf import read_vcf
import io
import os


def test_findcluster():
    """ Tests whether or not the proper clusters will be provided by the scipy clustering algorithm
    """

    cluster_dict = {
        'Reference': 3,
        'mysnpsSRR6683541': 1,
        'mysnpsSRR6683736': 1,
        'mysnpsSRR6683914': 1,
        'mysnpsSRR6683916': 2,
        'mysnpsSRR6683917': 4,
        'mysnpsSRR6703297': 3,
        'mysnpsSRR6967786': 1,
        'mysnpsSRR7128411': 3,
        'mysnpsSRR7128414': 3,
        'mysnpsSRR7130347': 3,
        'mysnpsSRR7130348': 3,
        'mysnpsSRR7130351': 3,
        'mysnpsSRR7130522': 3,
        'mysnpsSRR7135200': 5,
        'mysnpsSRR7167142': 3,
        'mysnpsSRR7168116': 3,
        'mysnpsSRR7168670': 3,
        'mysnpsSRR7170730': 3,
        'mysnpsSRR7221720': 3
    }

    data_frame = read_vcf("tests/data/core.vcf")

    result = find_clusters(data_frame)

    assert (result == cluster_dict)
