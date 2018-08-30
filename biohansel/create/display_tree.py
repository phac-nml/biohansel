import logging
import random
import os

from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

from typing import Dict

groups_dict_copy={}

def layout(node):
    if node.is_leaf():
        attr_face = AttrFace("name", fsize=50)
        
        faces.add_face_to_node(attr_face, node, 0, position="branch-right")
        group_name=faces.TextFace(groups_dict_copy[node.name], fsize=40)
        faces.add_face_to_node(group_name, node, column=0)
    
def get_colour():
    num1=random.randint(50,256)
    num2=random.randint(50,256)
    num3=random.randint(50,256)
    return "#{:02x}{:02x}{:02x}".format(num1, num2, num3)

def get_node_style():

    # Set dashed blue lines in all leaves
    node_style = NodeStyle()
    node_style["bgcolor"] = get_colour()
    
    return node_style

def display_tree(phylo_tree_path: str, groups_dict: Dict[str, str], output_folder_name:str) -> str:
    """ Displays the modified phylogenetic tree with subclades added to the original genome names, it also provides a colour-coded background based on subclade
    i.e. SRR238289-1.1.2
    Args:
        phylo_tree_path: the path to the user-defined phylogenetic tree
        groups_dict: the dictionary that contains the group information for each genome

    Returns:
        new_tree: the modified tree file to that displays the phylogenetic tree to the user
    """
   
    with open(phylo_tree_path) as file:
        logging.info(groups_dict)
        new_tree = Tree(file.read())
    
        unique_groups = list(set(groups_dict.values()))
        current_list = []

        for group in unique_groups:
            for genome, curr_group in groups_dict.items():
                if group == curr_group:
                    current_list.append(genome)
            curr_node=new_tree.get_common_ancestor(*current_list)
            curr_node.set_style(get_node_style())
            current_list=[]
        for genome, group in groups_dict.items():
            groups_dict_copy[genome]=group
        tree_style = TreeStyle()
        tree_style.layout_fn = layout
        tree_style.show_leaf_name = False
        tree_style.branch_vertical_margin=20
        tree_style.mode="c"
        tree_style.scale=120
        output_tree_path = os.path.join(output_folder_name, "phylo_tree.png")
        new_tree.render(output_tree_path, w=183, units="mm", tree_style=tree_style)
    return new_tree




