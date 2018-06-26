
import numpy as np
import argparse
from sys import argv
import math




def split(input_genomes):
    samples=[]
    with open(input_genomes) as my_file:
 
     for line in my_file:
        samples.append(line.strip())
    
    np.random.seed(42)
    indices=[]
    indices = np.random.permutation(samples)
    n_train = math.floor(len(indices)*0.75)
  
    
    test_indices=indices[n_train:]

    return test_indices


if __name__ == '__split__':
    split()