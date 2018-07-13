import math
from typing import List

import numpy as np


def  split_genomes(input_genomes: List[str])-> List[str]:
    """Extracts a set of genomes from the list of genomes inside the vcf file and creates a new list that contains
    the


    Args:
        input_genomes:

    Returns:
        test_indices:

    """
    samples = []
    for x in input_genomes:
        samples.append(x.strip())

    np.random.seed(42)
    indices = np.random.permutation(samples)
    n_train = math.floor(len(indices) * 0.75)

    test_indices = indices[n_train:]

    return test_indices
