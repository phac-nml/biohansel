import logging
import numpy as np
import pandas as pd
import scipy as sp

from typing import Dict, Set, List
from collections import defaultdict
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from biohansel.create.cluster import Cluster


def find_clusters(df: pd.DataFrame, group_size_range: tuple, distance_thresholds: List, pairwise_metric: str, linkage_method: str) -> Cluster:
    """
    Takes in a vcf file and creates clusters from the scipy hierarchy clustering algorithm
    
    Example:
    '/path/example.vcf' -> {'mysnpsSRR2323':'1', 'mysnpsSRR232323':'2'}
    
    Args:
        df: contains the vcf information in the DataFrame format
        min_group_size: the minimum child group size for each new subtype branching point from the parent group
        pairwise_metric: the distance metric used to calculate pairwise distances between SNVs
        linkage_method: the linkage method used to perform hierarchical clustering on SNVs

    Returns:
        cluster_dict: a dictionary indicating the cluster membership of each of the supplied genomes in
                the vcf file
    """
    


    distance_matrix = compute_distance_matrix(df, pairwise_metric )
    clustering_array = create_linkage_array(distance_matrix, linkage_method)

    flat_clusters = output_flat_clusters(clustering_array, df.columns, distance_matrix, group_size_range, distance_thresholds)
    cluster_result=Cluster(distance_matrix=distance_matrix, clustering_array= clustering_array, flat_clusters= flat_clusters)
    return cluster_result


def compute_distance_matrix(df: pd.DataFrame, pairwise_metric:str) -> np.ndarray:
    """Takes in a binary SNV state DataFrame and outputs a distance matrix using the hamming method

    Args:
        df: the DataFrame that contains only the samples' binary data

    Returns:
        a matrix of pair-wise distances between samples
    """

    return sp.spatial.distance.pdist(
        df.transpose(), metric=pairwise_metric)


def create_linkage_array(distance_matrix: np.ndarray, linkage_method: str) -> np.ndarray:
    """Takes in a distance matrix and outputs a hierarchical clustering linkage array

    Args:
        distance_matrix: a matrix of pair-wise distances between samples

    Returns:
        clustering_array: a hierarchical clustering linkage array that is calculated from the distance matrix
    """

    clustering_array = linkage(distance_matrix, method=linkage_method)

    return clustering_array


def output_flat_clusters(clustering_array: np.ndarray, genomes_only: List, distance_matrix: np.ndarray,
                         group_size_range: tuple, distance_thresholds: List) -> \
        pd.DataFrame:
    """Uses a set of thresholds to output a flat cluster of the linkage array

    Args:
        clustering_array: a hierarchical clustering linkage array that is calculated from the distance matrix
        genomes_only: an array of column names for the original SNV DataFrame
        distance_matrix: a matrix of pair-wise distances between samples
        group_size_range: the range of possivle child group sizes for each new subtype branching point from the parent group

    Returns:
         flat_clusters: an array of flat clusters from the clustering array
    """
    if distance_thresholds is not None:
        test_thresholds=distance_thresholds
    else:
        test_thresholds=np.sort(np.unique(distance_matrix))

    clusters = np.array([
        fcluster(clustering_array, t=distance, criterion='distance')
        for distance in test_thresholds])
    cluster_df = pd.DataFrame(np.array(clusters), index=test_thresholds)
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
            assign_hc_clusters(clusters_genomes_dict=clusters_genomes_dict, group_size_range=group_size_range))).fillna(
        '').loc[cluster_matrix.index, :]
    final_assigned_clusters = df_to_subtypes_dict(hierarchical_cluster_df)
    logging.debug(final_assigned_clusters)
    return final_assigned_clusters


def cluster_df_to_dict(df_clusters: pd.DataFrame) -> Dict[float, Dict[int, Set[str]]]:
    """Clusters dataframe to dict of threshold levels to clusters to members

    Args:
        df_clusters: The cluster dataframe with original cluster cluster_genomes

    Returns:
        clusters_genomes_dict: A dictionary of the genome cluster cluster_genomes with the thresholds used as the key
    """
    clusters_genomes_dict = {}
    for threshold in df_clusters.columns:
        clusters = df_clusters[threshold]  # type: pd.Series
        cluster_genomes = defaultdict(set)
        for genome, cluster in clusters.iteritems():
            cluster_genomes[cluster].add(genome)
        clusters_genomes_dict[threshold] = cluster_genomes

    return clusters_genomes_dict


def assign_hc_clusters(clusters_genomes_dict: Dict[float, Dict[int, Set[str]]], group_size_range: tuple) -> Dict[
    str, Dict[str, Set[str]]]:
    """
    Assigns the subtypes for each genome at each threshold level

    Args:
        clusters_genomes_dict: A dictionary of the genome cluster cluster_genomes with the thresholds used as the key
        group_size_range: the range of child group sizes for each new subtype branching point from the parent group

    Returns:
        output_clusters: A dictionary within a dictionary that contains each subgtype and the sets of genomes within
                        that subtype
    """
    min_group_size=group_size_range[0]
    max_group_size=group_size_range[1]
    output_subtypes = {threshold: {} for threshold in clusters_genomes_dict.keys()}
    sorted_thresholds = sorted(clusters_genomes_dict.keys())
    # initialize top level subtypes
    output_subtypes[sorted_thresholds[0]] = {str(cluster): genomes for cluster, genomes in
                                             clusters_genomes_dict[sorted_thresholds[0]].items()}

    for threshold_index in range(1, len(sorted_thresholds)):
        parent_hc_clusters = output_subtypes[sorted_thresholds[threshold_index - 1]]

        threshold = sorted_thresholds[threshold_index]
        cluster_genomes = clusters_genomes_dict[threshold]

        for subtype, parent_genomes in parent_hc_clusters.items():
            if len(parent_genomes) < min_group_size:
                continue
            subclade = 1
            for _, child_genomes in cluster_genomes.items():
                if len(child_genomes) < min_group_size or len(child_genomes)>max_group_size:
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
        cluster_dict: a dictionary that contains the cluster of each genome with each value being sets of
                    genomes

    Returns:
        out: A dictionary with values beings filled for all cells

    """

    out = {}
    for threshold, cluster_genomes in cluster_dict.items():
        out[threshold] = {genome: cluster for cluster, genomes in cluster_genomes.items() for genome in genomes}
    logging.debug(out)
    return out


def row_subtype(curr_genome: pd.Series) -> str:
    """
    Provides the last valid subtype for that particular genome, as that's the last valid branch point

    Args:
        curr_genome: the current subtypes for that particular genome

    Returns:
        max(row_unique.tolist(), key=len): the last valid subtype name for that genome

    """

    row_unique = curr_genome.unique()
    return max(row_unique.tolist(), key=len)
   


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
