from typing import List

import numpy as np
from scipy.cluster.hierarchy import to_tree, ClusterNode


def _scipy_tree_to_newick_list(node: ClusterNode, out: List[str], parent_distance: float, leaf_names: List[str]) -> List[str]:
    """Construct Newick tree from SciPy hierarchical clustering ClusterNode

    This is a recursive function to help build a Newick output string from a scipy.cluster.hierarchy.to_tree input with
    user specified leaf node names.

    Notes:
        This function is meant to be used with `to_newick`

    Args:
        node: Root node of tree from hierarchical clustering linkage matrix
        out: Newick string output accumulator list which needs to be reversed and concatenated (i.e. `''.join(out[::-1])`) to produce final Newick format string output
        parent_distance: Distance of parent node of `node`
        leaf_names: Leaf node names

    Returns:
        List of Newick output strings that need to be reversed and concatenated
    """
    if node.is_leaf():
        return out + [f'{leaf_names[node.id]}:{parent_distance - node.dist}']

    if len(out) > 0:
        out.append(f'):{parent_distance - node.dist}')
    else:
        out.append(');')
    out = _scipy_tree_to_newick_list(node.get_left(), out, node.dist, leaf_names)
    out.append(',')
    out = _scipy_tree_to_newick_list(node.get_right(), out, node.dist, leaf_names)
    out.append('(')
    return out


def cluster_node_to_newick(tree: ClusterNode, leaf_names: List[str]) -> str:
    """Newick tree output string from SciPy hierarchical clustering tree

    Convert a SciPy ClusterNode tree to a Newick format string.
    Use scipy.cluster.hierarchy.to_tree on a hierarchical clustering linkage matrix to create the root ClusterNode for the `tree` input of this function.

    Args:
        tree (scipy.cluster.hierarchy.ClusterNode): Output of scipy.cluster.hierarchy.to_tree from hierarchical clustering linkage matrix
        leaf_names (list of string): Leaf node names

    Returns:
        (string): Newick output string
    """
    newick_list = _scipy_tree_to_newick_list(tree, [], tree.dist, leaf_names)
    return ''.join(newick_list[::-1])


def newick_from_linkage(linkage_array: np.ndarray, leaf_names: List[str]) -> str:
    """Newick string from SciPy linkage array

    Hierarchical clustering of a distance matrix using SciPy's or fastcluster's linkage method returns a (N+1)x4 array.
    `scipy.cluster.hierarchy#to_tree` returns a `ClusterNode` object which can be traversed to generate a Newick format
    string.

    Args:
        linkage_array: Hierarchical clustering linkage array from SciPy's or `fastcluster`'s `linkage` method
        leaf_names: Tree leaf names (e.g. genome or sample names)

    Returns:
         Tree as Newick format string
    """
    tree_node = to_tree(linkage_array)
    return cluster_node_to_newick(tree_node, leaf_names)