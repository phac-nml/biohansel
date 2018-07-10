import scipy
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
import os
from ete3 import ClusterTree
import scipy.spatial.distance
import numpy as np
from itertools import combinations
from scipy import cluster
import pandas as pd
import io
from typing import Dict


def find_clusters(data_frame: pd.DataFrame) -> Dict[str, str]:
    """
    Takes in a vcf file and creates clusters from the scipy heirarchy clustering algorithm
    

    Example:
    '/path/example.vcf' -> {'mysnps2323':'1', 'mysnps232323':'2'}
    
    Args:
    data_frame: contains the vcf information in the DataFrame format

    Returns:
    clust_dict: a dictionary indicating the cluster membership of each of the supplied genomes in
                the vcf file
    data_frame: the data_frame version of the 

    """

    filtered_data_frame = data_frame.drop(['POS', 'REF', 'ALT'], 1)
    genomes_only = filtered_data_frame.columns
    distance_matrix = pdist(filtered_data_frame.transpose(), metric='hamming')

    clustering_array = linkage(distance_matrix, method='complete')
    thresholds = [0.2]

    clusters = np.array([
        fcluster(clustering_array, t=n, criterion='distance')
        for n in thresholds
    ])

    final = pd.DataFrame(
        np.array(clusters), index=thresholds, columns=genomes_only)

    subset = final.iloc[0]
    cluster_dict = subset.to_dict()

    return cluster_dict
