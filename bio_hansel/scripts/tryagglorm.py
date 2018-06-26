import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage

np.random.seed(123)
labels = ['ID_0','ID_1','ID_2','ID_3','ID_4']

x = np.corrcoef(np.random.random_sample([5,3])*10)
row_clusters = linkage(x, method='complete')    

def extract_levels(row_clusters, labels):
    clusters = {}
    for row in xrange(row_clusters.shape[0]):
        cluster_n = row + len(labels)
        # which clusters / labels are present in this row
        glob1, glob2 = row_clusters[row, 0], row_clusters[row, 1]

        # if this is a cluster, pull the cluster
        this_clust = []
        for glob in [glob1, glob2]:
            if glob > (len(labels)-1):
                this_clust += clusters[glob]
            # if it isn't, add the label to this cluster
            else:
                this_clust.append(glob)

        clusters[cluster_n] = this_clust
        print(clusters)
    return clusters

if __name__ == '__main__':
    extract_levels(row_clusters, labels)