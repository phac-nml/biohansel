import io
import os
import unittest
import warnings

import pytest
import pandas as pd
import numpy as np

from collections import defaultdict

from biohansel.create.io.parsers import parse_vcf
from biohansel.create.cluster_generator import find_clusters, compute_distance_matrix, create_linkage_array, \
    output_flat_clusters, cluster_df_to_dict, assign_hc_clusters, expand_sets, row_subtype, df_to_subtypes_dict

warnings.simplefilter(action='ignore', category=FutureWarning)


def test_findcluster():
    """ Tests whether or not the proper clusters will be provided by the scipy clustering algorithm
    """

    cluster_dict = {'Reference': '2.1.2',
                    'SRR6683541': '1',
                    'SRR6683736': '1',
                    'SRR6683914': '1',
                    'SRR6683916': '1',
                    'SRR6683917': '2.2',
                    'SRR6703297': '2.1',
                    'SRR6967786': '1',
                    'SRR7128411': '2.1.1.2',
                    'SRR7128414': '2.1.1.2',
                    'SRR7130347': '2.1.2',
                    'SRR7130348': '2.1.1.2',
                    'SRR7130351': '2.1.1.1',
                    'SRR7130522': '2.1.1.2',
                    'SRR7135200': '2.2',
                    'SRR7167142': '2.1.1.2',
                    'SRR7168116': '2.1',
                    'SRR7168670': '2.1.2',
                    'SRR7170730': '2.1.1.1',
                    'SRR7221720': '2.1.2'}

    _, data_frame = parse_vcf("tests/data/create/core.vcf")

    result = find_clusters(data_frame, 2)

    assert (result == cluster_dict)


def test_compute_distance_matrix():
    _, data_frame = parse_vcf("tests/data/create/test.vcf")
    test_array = compute_distance_matrix(data_frame)
    test_result = np.array([0.4, 0.4, 0.4, 0.8, 0, 0, 0.4, 0, 0.4, 0.4])

    np.testing.assert_array_equal(test_result, test_array)


def test_create_linkage_array():
    dataset = np.array([0.4, 0.4, 0.4, 0.8, 0, 0, 0.4, 0, 0.4, 0.4])
    expected_result = np.array(([1., 2., 0., 2.],
                                [3., 5., 0., 3., ],
                                [0., 6., 0.4, 4.],
                                [4., 7., 0.8, 5.]))
    test_result = create_linkage_array(dataset)
    np.testing.assert_array_equal(expected_result, test_result)


def test_output_flat_clusters():
    test_cluster_array = np.array(([1., 2., 0., 2.],
                                   [3., 5., 0., 3., ],
                                   [0., 6., 0.4, 4.],
                                   [4., 7., 0.8, 5.]))
    test_distance_array = np.array([0.4, 0.4, 0.4, 0.8, 0, 0, 0.4, 0, 0.4, 0.4])
    test_genomes = ['Reference', 'SRR6683541', 'SRR6683736', 'SRR6683914', 'SRR6683916']
    group_size = 2

    test_result = output_flat_clusters(test_cluster_array, test_genomes, test_distance_array, group_size)

    expected_clusters = {'Reference': '1',
                         'SRR6683541': '1',
                         'SRR6683736': '1',
                         'SRR6683914': '1',
                         'SRR6683916': '2'}

    assert (test_result == expected_clusters)


def test_cluster_df_to_dict():
    test_cluster_matrix = pd.read_csv("tests/data/create/test_dataframe.csv", index_col=0)
    test_dict = {1: defaultdict(set,
                                {1: {'SRR6683541',
                                     'SRR6683736',
                                     'SRR6683914',
                                     'SRR6683916',
                                     'SRR6967786'},
                                 2: {'Reference',
                                     'SRR6683917',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7135200',
                                     'SRR7167142',
                                     'SRR7168116',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'}}),
                 2: defaultdict(set,
                                {1: {'SRR6683541',
                                     'SRR6683736',
                                     'SRR6683914',
                                     'SRR6683916',
                                     'SRR6967786'},
                                 2: {'Reference',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168116',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 3: {'SRR6683917', 'SRR7135200'}}),
                 3: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736', 'SRR6683914', 'SRR6967786'},
                                 2: {'SRR6683916'},
                                 3: {'Reference',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168116',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 4: {'SRR6683917', 'SRR7135200'}}),
                 4: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736', 'SRR6683914', 'SRR6967786'},
                                 2: {'SRR6683916'},
                                 3: {'Reference',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168116',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 4: {'SRR6683917'},
                                 5: {'SRR7135200'}}),
                 5: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736', 'SRR6967786'},
                                 2: {'SRR6683914'},
                                 3: {'SRR6683916'},
                                 4: {'Reference',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168116',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 5: {'SRR6683917'},
                                 6: {'SRR7135200'}}),
                 6: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736', 'SRR6967786'},
                                 2: {'SRR6683914'},
                                 3: {'SRR6683916'},
                                 4: {'Reference',
                                     'SRR6703297',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 5: {'SRR7168116'},
                                 6: {'SRR6683917'},
                                 7: {'SRR7135200'}}),
                 7: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736', 'SRR6967786'},
                                 2: {'SRR6683914'},
                                 3: {'SRR6683916'},
                                 4: {'Reference',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 5: {'SRR6703297'},
                                 6: {'SRR7168116'},
                                 7: {'SRR6683917'},
                                 8: {'SRR7135200'}}),
                 8: defaultdict(set,
                                {1: {'SRR6683541', 'SRR6683736'},
                                 2: {'SRR6967786'},
                                 3: {'SRR6683914'},
                                 4: {'SRR6683916'},
                                 5: {'Reference',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 6: {'SRR6703297'},
                                 7: {'SRR7168116'},
                                 8: {'SRR6683917'},
                                 9: {'SRR7135200'}}),
                 9: defaultdict(set,
                                {1: {'SRR6683541'},
                                 2: {'SRR6683736'},
                                 3: {'SRR6967786'},
                                 4: {'SRR6683914'},
                                 5: {'SRR6683916'},
                                 6: {'Reference',
                                     'SRR7128411',
                                     'SRR7128414',
                                     'SRR7130347',
                                     'SRR7130348',
                                     'SRR7130351',
                                     'SRR7130522',
                                     'SRR7167142',
                                     'SRR7168670',
                                     'SRR7170730',
                                     'SRR7221720'},
                                 7: {'SRR6703297'},
                                 8: {'SRR7168116'},
                                 9: {'SRR6683917'},
                                 10: {'SRR7135200'}}),
                 10: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7128411',
                                      'SRR7128414',
                                      'SRR7130348',
                                      'SRR7130351',
                                      'SRR7130522',
                                      'SRR7167142',
                                      'SRR7170730'},
                                  7: {'Reference', 'SRR7130347', 'SRR7168670', 'SRR7221720'},
                                  8: {'SRR6703297'},
                                  9: {'SRR7168116'},
                                  10: {'SRR6683917'},
                                  11: {'SRR7135200'}}),
                 11: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7128411',
                                      'SRR7128414',
                                      'SRR7130348',
                                      'SRR7130351',
                                      'SRR7130522',
                                      'SRR7167142',
                                      'SRR7170730'},
                                  7: {'SRR7130347', 'SRR7168670', 'SRR7221720'},
                                  8: {'Reference'},
                                  9: {'SRR6703297'},
                                  10: {'SRR7168116'},
                                  11: {'SRR6683917'},
                                  12: {'SRR7135200'}}),
                 12: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7128411',
                                      'SRR7128414',
                                      'SRR7130348',
                                      'SRR7130351',
                                      'SRR7130522',
                                      'SRR7167142',
                                      'SRR7170730'},
                                  7: {'SRR7130347', 'SRR7168670'},
                                  8: {'SRR7221720'},
                                  9: {'Reference'},
                                  10: {'SRR6703297'},
                                  11: {'SRR7168116'},
                                  12: {'SRR6683917'},
                                  13: {'SRR7135200'}}),
                 13: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7130351', 'SRR7170730'},
                                  7: {'SRR7128411',
                                      'SRR7128414',
                                      'SRR7130348',
                                      'SRR7130522',
                                      'SRR7167142'},
                                  8: {'SRR7130347', 'SRR7168670'},
                                  9: {'SRR7221720'},
                                  10: {'Reference'},
                                  11: {'SRR6703297'},
                                  12: {'SRR7168116'},
                                  13: {'SRR6683917'},
                                  14: {'SRR7135200'}}),
                 14: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7130351', 'SRR7170730'},
                                  7: {'SRR7128411', 'SRR7128414', 'SRR7130348', 'SRR7167142'},
                                  8: {'SRR7130522'},
                                  9: {'SRR7130347'},
                                  10: {'SRR7168670'},
                                  11: {'SRR7221720'},
                                  12: {'Reference'},
                                  13: {'SRR6703297'},
                                  14: {'SRR7168116'},
                                  15: {'SRR6683917'},
                                  16: {'SRR7135200'}}),
                 15: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7130351'},
                                  7: {'SRR7170730'},
                                  8: {'SRR7128411', 'SRR7128414', 'SRR7167142'},
                                  9: {'SRR7130348'},
                                  10: {'SRR7130522'},
                                  11: {'SRR7130347'},
                                  12: {'SRR7168670'},
                                  13: {'SRR7221720'},
                                  14: {'Reference'},
                                  15: {'SRR6703297'},
                                  16: {'SRR7168116'},
                                  17: {'SRR6683917'},
                                  18: {'SRR7135200'}}),
                 16: defaultdict(set,
                                 {1: {'SRR6683541'},
                                  2: {'SRR6683736'},
                                  3: {'SRR6967786'},
                                  4: {'SRR6683914'},
                                  5: {'SRR6683916'},
                                  6: {'SRR7130351'},
                                  7: {'SRR7170730'},
                                  8: {'SRR7128414', 'SRR7167142'},
                                  9: {'SRR7128411'},
                                  10: {'SRR7130348'},
                                  11: {'SRR7130522'},
                                  12: {'SRR7130347'},
                                  13: {'SRR7168670'},
                                  14: {'SRR7221720'},
                                  15: {'Reference'},
                                  16: {'SRR6703297'},
                                  17: {'SRR7168116'},
                                  18: {'SRR6683917'},
                                  19: {'SRR7135200'}})}
    test_result = cluster_df_to_dict(test_cluster_matrix)

    assert (len(test_dict) == len(test_result))


def test_assign_hc_clusters():
    test_dict = {1: defaultdict(set, {1: {'SRR6683914', 'Reference', 'SRR6683541', 'SRR6683736'}, 2: {'SRR6683916'}}),
                 2: defaultdict(set,
                                {1: {'SRR6683914', 'SRR6683541', 'SRR6683736'}, 2: {'Reference'}, 3: {'SRR6683916'}})}
    group_size = 2
    expected_result = {1: {'1': {'SRR6683914', 'Reference', 'SRR6683541', 'SRR6683736'}, '2': {'SRR6683916'}},
                       2: {'1': {'SRR6683914', 'SRR6683541', 'SRR6683736'}}}
    test_result = assign_hc_clusters(test_dict, group_size)

    assert (test_result == expected_result)


def test_expand_sets():
    test_dict = {1: {'1': {'SRR6683914', 'Reference', 'SRR6683541', 'SRR6683736'}, '2': {'SRR6683916'}},
                 2: {'1': {'SRR6683914', 'SRR6683541', 'SRR6683736'}}}
    expected_result = {
        1: {'Reference': '1', 'SRR6683914': '1', 'SRR6683541': '1', 'SRR6683736': '1', 'SRR6683916': '2'},
        2: {'SRR6683914': '1', 'SRR6683541': '1', 'SRR6683736': '1'}}

    test_result = expand_sets(test_dict)

    assert (test_result == expected_result)


def test_row_subtype():
    test_data = pd.Series(data=[1, 1, 2])
    expected_result = 2
    test_result = row_subtype(test_data)

    assert (test_result == expected_result)


def test_df_to_subtypes_dict():
    test_data = pd.read_csv("tests/data/create/cluster_genomes_test.csv", index_col=0)
    expected_result = {'SRR6683541': 1, 'SRR6683736': 1, 'SRR6683914': 1, 'Reference': 1, 'SRR6683916': 2}
    test_result = df_to_subtypes_dict(test_data)

    assert (expected_result == test_result)
