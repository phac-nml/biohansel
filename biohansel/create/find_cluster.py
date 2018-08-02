from typing import Dict, Set, List, Union

import logging
import numpy as np
import pandas as pd
import scipy as sp
import pprint

from collections import defaultdict
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


def find_clusters(df: pd.DataFrame, min_group_size: int) -> Dict[str, int]:
    """
    Takes in a vcf file and creates clusters from the scipy hierarchy clustering algorithm
    
    Example:
    '/path/example.vcf' -> {'mysnpsSRR2323':'1', 'mysnpsSRR232323':'2'}
    
    Args:
        df: contains the vcf information in the DataFrame format
        min_group_size: the minimum child group size for each new subtype branching point from the parent group

    Returns:
        cluster_dict: a dictionary indicating the cluster membership of each of the supplied genomes in
                the vcf file
    """

    distance_matrix = compute_distance_matrix(df)
    clustering_array = create_linkage_array(distance_matrix)

    flat_clusters = output_flat_clusters(clustering_array, df.columns, distance_matrix, min_group_size)

    return flat_clusters


def compute_distance_matrix(df: pd.DataFrame) -> List:
    """Takes in a binary SNV state DataFrame and outputs a distance matrix using the hamming method

    Args:
        df: the DataFrame that contains only the samples' binary data

    Returns:
        a matrix of pair-wise distances between samples
    """
    return sp.spatial.distance.pdist(
        df.transpose(), metric='hamming')


def create_linkage_array(distance_matrix: List) -> List:
    """Takes in a distance matrix and outputs a hierarchical clustering linkage array

    Args:
        distance_matrix: a matrix of pair-wise distances between samples

    Returns:
        clustering_array: a hierarchical clustering linkage array that is calculated from the distance matrix
    """
    clustering_array = linkage(distance_matrix, method='complete')
    return clustering_array


def output_flat_clusters(clustering_array: List, genomes_only: List, distance_matrix: List, min_group_size: int) -> \
        pd.DataFrame:
    """Uses a set of thresholds to output a flat cluster of the linkage array

    Args:
        clustering_array: a hierarchical clustering linkage array that is calculated from the distance matrix
        genomes_only: an array of column names for the original SNV DataFrame
        distance_matrix: a matrix of pair-wise distances between samples

    Returns:
         flat_clusters: an array of flat clusters from the clustering array
    """
    clusters = np.array([
        fcluster(clustering_array, t=distance, criterion='distance')
        for distance in np.sort(np.unique(distance_matrix))])
    cluster_df = pd.DataFrame(np.array(clusters), index=np.sort(np.unique(distance_matrix)))
    matrix = cluster_df.transpose()
    matrix.index = genomes_only
    matrix = matrix.sort_values(
        by=[cluster_grouping for cluster_grouping in matrix.columns.sort_values(ascending=False)])
    thresholds = [threshold for threshold in matrix.columns.sort_values(ascending=False)]

    cluster_dict = defaultdict(list)

    for threshold in thresholds:
        grouping = '-'.join([str(cluster) for cluster in matrix[threshold]])
        cluster_dict[grouping].append(threshold)

    cluster_matrix = matrix[[thresholds[0] for thresholds in cluster_dict.values()]]
    cluster_matrix.columns = [index for index, thresholds in enumerate(cluster_matrix.columns)]

    cluster_matrix = cluster_matrix.drop([0], axis=1)

    clusters_genomes_dict = cluster_df_to_dict(cluster_matrix)
    hierarchical_cluster_df = pd.DataFrame(
        expand_sets(
            assign_hc_clusters(clusters_genomes_dict=clusters_genomes_dict, min_group_size=min_group_size))).fillna(
        '').loc[cluster_matrix.index, :]
    final_assigned_clusters = df_to_subtypes_dict(hierarchical_cluster_df)

    logging.debug(pprint.pformat(final_assigned_clusters))

    return final_assigned_clusters


def cluster_df_to_dict(df_clusters: pd.DataFrame) -> Dict[float, Dict[int, Set[str]]]:
    """Clusters dataframe to dict of threshold levels to clusters to members

    Args:
        df_clusters: The cluster dataframe with original cluster groupings

    Returns:
        clusters_genomes_dict: A dictionary of the genome cluster groupings with the thresholds used as the key
    """
    clusters_genomes_dict = {}
    for threshold in df_clusters.columns:
        clusters = df_clusters[threshold]  # type: pd.Series
        cluster_genomes = defaultdict(set)
        for genome, cluster in clusters.iteritems():
            cluster_genomes[cluster].add(genome)
        clusters_genomes_dict[threshold] = cluster_genomes
    return clusters_genomes_dict


def assign_hc_clusters(clusters_genomes_dict: Dict[float, Dict[int, Set[str]]], min_group_size: int) -> Dict[
    str, Dict[str, Set[str]]]:
    """
    Assigns the subtypes for each genome at each threshold level

    Args:
        clusters_genomes_dict: A dictionary of the genome cluster groupings with the thresholds used as the key
        min_group_size: the minimum child group size for each new subtype branching point from the parent group

    Returns:
        output_clusters: A dictionary within a dictionary that contains each subgtype and the sets of genomes within
                        that subtype
    """

    output_subtypes = {threshold: {} for threshold in clusters_genomes_dict.keys()}
    sorted_thresholds = sorted(clusters_genomes_dict.keys())
    # initialize top level subtypes
    output_subtypes[sorted_thresholds[0]] = {str(cluster): genomes for cluster, genomes in clusters_genomes_dict[sorted_thresholds[0]].items()}

    for threshold_index in range(1, len(sorted_thresholds)):
        parent_hc_clusters = output_subtypes[sorted_thresholds[threshold_index - 1]]

        threshold = sorted_thresholds[threshold_index]
        cluster_genomes = clusters_genomes_dict[threshold]

        for subtype, parent_genomes in parent_hc_clusters.items():
            if len(parent_genomes) < min_group_size:
                continue
            subclade = 1
            for _, child_genomes in cluster_genomes.items():
                if len(child_genomes) < min_group_size:
                    continue

                if parent_genomes == child_genomes:
                    output_subtypes[threshold][subtype] = parent_genomes
                elif child_genomes.issubset(parent_genomes):
                    new_subgroup_name = f'{subtype}.{subclade}'
                    if len(parent_genomes - child_genomes) >= min_group_size:
                        output_subtypes[threshold][new_subgroup_name] = child_genomes
                    else:
                        output_subtypes[threshold][subtype] = child_genomes
                    subclade += 1
    return output_subtypes


def expand_sets(cluster_dict: Dict[str, Dict[str, Set[str]]]) -> Dict[str, Dict[str, str]]:
    """
    Fills in the rest of the cells in the dataframe that had previously been unfilled

    Args:
        cluster_dict: a dictionary that contains the cluster groupings of each genome with each value being sets of
                    genomes

    Returns:
        modified_cluster_dict: A dictionary with values beings filled for all cells

    """
    modified_cluster_dict = {}
    for threshold, groupings in cluster_dict.items():
        threshold_dict = {}
        for grouping, genomes in groupings.items():
            for genome in genomes:
                threshold_dict[genome] = grouping
        modified_cluster_dict[threshold] = threshold_dict
    return modified_cluster_dict


def row_subtype(curr_genome: pd.Series) -> str:
    """
    Provides the last valid subtype for that particular genome, as that's the last valid branch point

    Args:
        curr_genome: the current subtypes for that particular genome

    Returns:
        row_unique[-1]: the last valid subtype name for that genome

    """
    row_unique = curr_genome.unique()
    row_unique = row_unique[(row_unique != '') & (~pd.isnull(row_unique))]
    return row_unique[-1]


def df_to_subtypes_dict(cluster_genomes_df: pd.DataFrame) -> Dict[str, str]:
    """
    Provides a dictionary of the subtypes with the genome being the key and the subtype being the value

    Args:
        cluster_genomes_df: the dataframe containing the genomes and their respective subtypes at each threshold

    Returns:
        final_cluster_dict: the dictionary that contains the final subtype assignments for each genome

    """
    final_cluster_dict = {}
    for genome, row in cluster_genomes_df.apply(row_subtype, axis=1).iteritems():
        final_cluster_dict[genome] = row
    return final_cluster_dict
