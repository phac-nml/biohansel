import logging

from ete3 import Tree

from typing import Dict




def display_tree(phylo_tree_path: str, groups_dict: Dict[str, str]) -> str:
    """ Displays the modified phylogenetic tree with subclades added to the original genome names
    i.e. SRR238289-1.1.2
    Args:
        phylo_tree_path: the path to the user-defined phylogenetic tree
        groups_dict: the dictionary that contains the group information for each genome

    Returns:
        new_tree: the modified tree file to that displays the phylogenetic tree to the user
    """
   
    with open(phylo_tree_path) as file:
        new_tree = file.read()
        for genome, group in groups_dict.items():
            new_name = f"{genome}-{group}"
            new_tree = new_tree.replace(genome, new_name)
        tree_diagram = Tree(new_tree)
        logging.info(f"{tree_diagram}\n")
    
    return new_tree
