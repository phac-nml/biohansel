from Bio import Phylo
import os

home_folder=os.path.expanduser('~')
tree=Phylo.read(f"{home_folder}/ecoli-genomes/e-coliO157-tree_file.tree",'newick')
print (tree)