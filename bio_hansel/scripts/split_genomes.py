import numpy as np
import argparse
from sys import argv
import math


def split_genomes(input_genomes: list):
    samples = []
    for x in input_genomes:
        samples.append(x.strip())

    np.random.seed(42)
    indices = []
    indices = np.random.permutation(samples)
    n_train = math.floor(len(indices) * 0.75)

    test_indices = indices[n_train:]

    return test_indices
