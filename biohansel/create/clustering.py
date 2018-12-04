import logging
from typing import Union, List

import attr
import numpy as np
import pandas as pd
from fastcluster import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import pdist

from biohansel.create.subtype_node import SubtypeNode


@attr.s
class HClust(object):
    """Class to store Hierarchical Clustering and flat cluster information on some input matrix.

    Attributes:
        distance_matrix: Distance matrix computed from input matrix
        pdist_metric: Pairwise distance metric used to compute distance matrix (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
        linkage_array: Hierarchical clustering linkage array computed by `fastcluster` `linkage`
        linkage_method: Hierarchical clustering linkage method used to produce `linkage_array` (see http://danifold.net/fastcluster.html?section=3)
        df_clusters: Table of flat clusters defined at multiple levels from greatest to least distance threshold.
        distance_thresholds: Distance thresholds used to define flat clusters.
    """
    df_binary_snvs: pd.DataFrame = attr.ib(default=None)
    distance_matrix: np.ndarray = attr.ib(default=None)
    pdist_metric: str = attr.ib(default='euclidean')
    linkage_array: np.ndarray = attr.ib(default=None)
    linkage_method: str = attr.ib(default='single')
    df_clusters: pd.DataFrame = attr.ib(default=None)
    distance_thresholds: Union[List[float], np.ndarray] = attr.ib(default=None)

    @classmethod
    def from_binary_snp_matrix(cls,
                               df_binary_snvs: pd.DataFrame,
                               distance_thresholds: List[float] = None,
                               pdist_metric: str = 'euclidean',
                               linkage_method: str = 'single') -> 'HClust':
        """Hierarchical clustering and flat clusters from a binary matrix of SNP absence/presence.

        Args:
            df_binary_snvs: Binary SNV absence/presence matrix
            distance_thresholds: Optional distance thresholds to define flat clusters at
            pdist_metric: Pairwise distance metric to use (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
            linkage_method: Hierarchical clustering linkage method to use (see http://danifold.net/fastcluster.html?section=3)

        Returns:

        """
        obj = cls()
        obj.df_binary_snvs = df_binary_snvs
        obj.pdist_metric = pdist_metric
        obj.linkage_method = linkage_method

        logging.info(f'Computing distance matrix using the "{pdist_metric}" metric from binary SNV matrix of shape '
                     f'{df_binary_snvs.shape}.')
        obj.distance_matrix = pdist(df_binary_snvs.transpose().values, metric=pdist_metric)
        logging.info(f'Computed distance matrix of size {obj.distance_matrix.size} using the "{pdist_metric}" metric '
                     f'from binary SNV matrix of shape {df_binary_snvs.shape}')
        logging.info(f'Performing hierarchical clustering of distance matrix using "{linkage_method}" linkage method.')
        obj.linkage_array = linkage(obj.distance_matrix, method=linkage_method)
        logging.info(f'Performed hierarchical clustering using "{linkage_method}" linkage method to produce linkage '
                     f'array with dimensions {obj.linkage_array.shape}.')
        obj._set_distance_thresholds(distance_thresholds)
        obj.df_clusters = fclusters_at_multiple_distances(obj.linkage_array, obj.distance_thresholds,
                                                          df_binary_snvs.columns)
        logging.info(f'Defined flat clusters at N={obj.df_clusters.shape[1]} non-redundant distance threshold levels.')

        return obj

    def subtypes_from_clusters(self, min_group_size: int = 5, min_num_snvs: int = 2, start_subtype: int = 1):
        """Create subtype membership tree from clusters table with minimum subtype membership size constraint.

        Args:
            min_group_size: Minimum subtype membership size
            min_num_snvs: Minimum number of SNVs to differentiate a subgroup from everything else
            start_subtype: Subtype designation to start at. If `0`, then the first subtype will be `0` and the next `1`, `2`, and so on.

        Returns:
            Root `SubtypeNode` for tree of hierachical subtypes and their info.
        """
        return SubtypeNode.from_fclusters(df_clusters=self.df_clusters,
                                          min_group_size=min_group_size,
                                          start_subtype=start_subtype,
                                          min_num_snvs=min_num_snvs,
                                          df_bin=self.df_binary_snvs)

    def _set_distance_thresholds(self, distance_thresholds: Union[List[float], np.ndarray] = None) -> None:
        if distance_thresholds:
            self.distance_thresholds = distance_thresholds
            logging.info(f'n={len(self.distance_thresholds)} distances specified')
        else:
            self.distance_thresholds = np.sort(np.unique(self.distance_matrix))
            logging.info(f'No distance thresholds specified. Using n={self.distance_thresholds.size} unique distances '
                         f'from the distance matrix')


def fclusters_at_multiple_distances(linkage_array: np.ndarray,
                                    distance_thresholds: List[float],
                                    index_names: Union[List, pd.Series]) -> pd.DataFrame:
    """Flat clusters defined at multiple distance levels from hierarchical clustering linkage array.

    Only unique cluster profiles are returned from this method so if the clusters produced at distance 0.1 and distance
    0.2 are the same, only the clusters at 0.1

    Args:
        linkage_array: Hierarchical clustering linkage array from `fastcluster`/SciPy `linkage`
        distance_thresholds: Distance thresholds to define flat clusters at
        index_names: Index names (e.g. genome names)

    Returns:
        Table of non-redundant clusters defined at multiple levels sorted by greatest distance to least distance threshold
    """
    # sort by least to greatest distance
    distance_thresholds = np.sort(distance_thresholds)[::-1]
    unique_clusters = set()
    threshold_cluster = {}
    for threshold in distance_thresholds:
        clusters = fcluster(linkage_array, t=threshold, criterion='distance')
        n_unique_clusters = np.unique(clusters).size
        if n_unique_clusters == 1:
            continue
        clusters_string = clusters.tostring()
        # unique clusters produced at this threshold?
        if clusters_string not in unique_clusters:
            unique_clusters.add(clusters_string)
            threshold_cluster[threshold] = clusters
    df = pd.DataFrame(threshold_cluster, index=index_names)
    # sort columns by greatest distance threshold to least
    df = df[np.sort(df.columns)[::-1]]
    return df.sort_values(by=df.columns.tolist(), ascending=True)
