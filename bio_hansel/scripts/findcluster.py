import scipy
from scipy.cluster.hierarchy import fcluster
import os
from ete3 import ClusterTree
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance
import numpy as np
from itertools import combinations
from scipy import cluster

def findClusters(tree_file, input_genomes, output_directory):
    """

    Returns a groups file that indicates the group membership of each of the genomes based on their tree distance
    """
    
    # tree_file=f"{homedirectory}/ecoli-genomes/e-coliO157-tree_file.tree"
    # input_genomes=f"{homedirectory}/ecoli-genomes/List_EcoliO157_genomes.txt"
    with open(tree_file, 'r') as file:
        with open(input_genomes,'r') as input_file:
            tree_input=file.read().replace('\n', '')
            tree = ClusterTree(tree_input)
            leaves = tree.get_leaf_names()
            idx_dict={}
            for i in range(len(leaves)):
                idx_dict[leaves[i]]=i
                
        # idx_dict = {'A':0,'B':1,'C':2,'D':3, 'E':4}
        # values=list(idx_dict.values())
        # keys=list(idx_dict.keys())
        # idx_labels = [keys[values[i]] for i in range(0, len(idx_dict))]

        dmat = np.zeros((len(leaves),len(leaves)))

        for l1,l2 in combinations(leaves,2):
            d = tree.get_distance(l1,l2)
            dmat[idx_dict[l1],idx_dict[l2]] = dmat[idx_dict[l2],idx_dict[l1]] = d

        # print ('Distance:')
        # print (dmat)


        schlink = sch.linkage(scipy.spatial.distance.squareform(dmat),method='average',metric='euclidean')
        clusters=fcluster(schlink,t=0.2, criterion='distance')

        # print ('Linkage from scipy:')
        # print (schlink)

        print ('Clusters:')
        cluster_dict={}

        for i in range(len(leaves)):
                cluster_dict[leaves[i]]=clusters[i]
        
        with open('newfile.txt', 'w+') as groupsfile:
            groupsfile.write("genomes\tgroup\n")
            for key, value in cluster_dict.items():
                groupsfile.write(f"{key}\t{value}\n")
                
        
        print (cluster_dict)
        return groupsfile
findClusters(3,2,3)