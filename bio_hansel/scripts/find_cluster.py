from typing import List, Dict

import logging
import numpy as np
import pandas as pd
import scipy as sp

from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


def find_clusters(df: pd.DataFrame, min_threshold: float,
                  max_threshold: float) -> Dict[str, int]:
    """
    Takes in a vcf file and creates clusters from the scipy hierarchy clustering algorithm
    
    Example:
    '/path/example.vcf' -> {'mysnpsSRR2323':'1', 'mysnpsSRR232323':'2'}
    
    Args:
        df: contains the vcf information in the DataFrame format
        min_threshold: minimum threshold for the flat cluster array
        max_threshold: maximum threshold for the flat cluster array

    Returns:
        cluster_dict: a dictionary indicating the cluster membership of each of the supplied genomes in
                the vcf file
    """

    distance_matrix = compute_distance_matrix(df)
    clustering_array = create_linkage_array(distance_matrix)

    flat_clusters = output_flat_clusters(clustering_array, df.columns,
                                         min_threshold, max_threshold)

    subset = flat_clusters.iloc[0]
    cluster_dict = subset.to_dict()

    logging.debug(cluster_dict)
    return cluster_dict


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


def output_flat_clusters(clustering_array: List, genomes_only: List, min_threshold: float, max_threshold: float) -> \
        pd.DataFrame:
    """Uses a set of thresholds to output a flat cluster of the linkage array

    Args:
        clustering_array: a hierarchical clustering linkage array that is calculated from the distance matrix
        genomes_only: an array of column names for the original SNV DataFrame
        min_threshold: the minimum threshold for the fcluster array
        max_threshold: the maximum threshold for the fcluster array


    Returns:
         flat_clusters: an array of flat clusters from the clustering array
    """

    thresholds = []
    max_threshold = 1
    min_threshold = 0.3
    numbers = int((max_threshold - min_threshold) / 0.2)
    numbers += 1
    curr_number = min_threshold
    for x in range(0, numbers):
        thresholds.append(curr_number)
        curr_number = curr_number + 0.2

    clusters = np.array([
        fcluster(clustering_array, t=n, criterion='distance')
        for n in thresholds
    ])

    flat_clusters = pd.DataFrame(
        np.array(clusters), index=thresholds, columns=genomes_only)

    return flat_clusters
