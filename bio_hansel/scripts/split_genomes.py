
import numpy as np
import argparse
from sys import argv
import math




def main():
    print('hello world')
    samples=[]
    with open('List_EcoliO157_genomes.txt') as my_file:
     for line in my_file:
        samples.append(line.strip())
    
    print(samples)

   
    np.random.seed(42)
    indices=[]
    indices = np.random.permutation(samples)
    n_train = math.floor(len(indices)*0.75)
    print(n_train)
    train_indices, test_indices = indices[:n_train], indices[n_train:]
    print(train_indices)
    print(test_indices)


if __name__ == '__main__':
    main()