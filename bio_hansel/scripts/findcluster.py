import scipy
from scipy.cluster.hierarchy import fcluster
import os
from ete3 import ClusterTree
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance
import numpy as np
from itertools import combinations
from scipy import cluster

def findClusters(tree_file: str, input_genomes: str, output_directory: str):
    """

    Returns a groups file that indicates the group membership of each of the genomes based on their tree distance
    """
   
    with open(tree_file, 'r') as file:
        with open(input_genomes,'r') as input_file:
            tree_input=file.read().replace('\n', '')
            tree = ClusterTree(tree_input)
            leaves = tree.get_leaf_names()
            idx_dict={}
            for i in range(len(leaves)):
                idx_dict[leaves[i]]=i
                
      
        dmat = np.zeros((len(leaves),len(leaves)))

        for l1,l2 in combinations(leaves,2):
            d = tree.get_distance(l1,l2)
            dmat[idx_dict[l1],idx_dict[l2]] = dmat[idx_dict[l2],idx_dict[l1]] = d

        schlink = sch.linkage(scipy.spatial.distance.squareform(dmat),method='average',metric='euclidean')
        clusters=fcluster(schlink,t=0.2, criterion='distance')

       

        print ('Clusters:')
        cluster_dict={}

        for i in range(len(leaves)):
                cluster_dict[leaves[i]]=clusters[i]
        
        cluster_filepath=f"{output_directory}/clusterfile.txt"
        
        with open(cluster_filepath, 'w+') as clusterfile:
            clusterfile.write("genomes\tgroup\n")
            for key, value in cluster_dict.items():
                clusterfile.write(f"{key}\t{value}\n")
                
        
        print (cluster_dict)
        return cluster_filepath
