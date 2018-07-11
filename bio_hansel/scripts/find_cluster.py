from typing import Dict

import numpy as np
import pandas as pd
import scipy as sp

from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist



def find_clusters(df: pd.DataFrame) -> Dict[str, str]:
    """
    Takes in a vcf file and creates clusters from the scipy hierarchy clustering algorithm
    
    Example:
    '/path/example.vcf' -> {'mysnps2323':'1', 'mysnps232323':'2'}
    
    Args:
    df: contains the vcf information in the DataFrame format

    Returns:
    cluster_dict: a dictionary indicating the cluster membership of each of the supplied genomes in
                the vcf file

    """

    filtered_df = filter_df(df)
    distance_matrix = compute_distance_matrix(filtered_df)
    clustering_array = create_linkage_array(distance_matrix)

    flat_clusters = output_flat_clusters(clustering_array, filtered_df.columns)

    subset = flat_clusters.iloc[0]
    cluster_dict = subset.to_dict()
    return cluster_dict

#Still need to add docstrings and to also add the variable types to each individual function
def filter_df(df: pd.DataFrame):
    filtered_df = df.drop(['POS', 'REF', 'ALT'], 1)
    return filtered_df


def compute_distance_matrix(filtered_df: pd.DataFrame):
    distance_matrix = sp.spatial.distance.pdist(filtered_df.transpose(), metric='hamming')
    return distance_matrix


def create_linkage_array(distance_matrix):
    clustering_array = linkage(distance_matrix, method='complete')
    return clustering_array


def output_flat_clusters(clustering_array:list, genomes_only):
    thresholds = [0.2, 0.3, 0.6]

    clusters = np.array([
        fcluster(clustering_array, t=n, criterion='distance')
        for n in thresholds
    ])

    flat_clusters = pd.DataFrame(
        np.array(clusters), index=thresholds, columns=genomes_only)

    return flat_clusters
