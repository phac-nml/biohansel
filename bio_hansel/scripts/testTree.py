from ete3 import Tree, PhyloTree, ClusterTree

import os

home_folder=os.path.expanduser('~')

with open(f"{home_folder}/ecoli-genomes/e-coliO157-tree_file.tree", 'r') as myfile:
    data=myfile.read().replace('\n', '')

t=ClusterTree(data)


# t.set_species_naming_function(lambda node: node.name.split("_")[0] )

# print(t.get_ascii(attributes=["name", "species"], show_internal=False ))

# for node in t.split_by_dups():
#     print (node)
print(t)