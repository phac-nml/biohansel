from scipy.signal import argrelextrema, savgol_filter
from typing import List
import numpy as np


'''
[calc_avg_kmer_depth]
    Input: A dictionary containing the result of running Jellyfish's histo method in the k-mer counts
    Output: The estimated average kmer coverage value
    Desc: To find the average kmer coverage value, we use the assumption that most data follows the following curve:
    |.
    | .       |------- Here is the average k-mer coverage depth value.
    | .     . | .
    | .    .  |  .
    | .   .   |   .
    |  . .    |    . .
    |   .     |       . . .
    |______________________________
        
    Where the maxima's value would be your average k-mer coverage value by accessing it's key. 
'''


def calc_avg_kmer_depth(kmers: dict) -> int:
    maximum = find_maxima(apply_savgol_filt(list(kmers.values())))

    for key, value in kmers.items():
        if maximum == value:
            kmer_coverage = key

    return kmer_coverage


def calc_error_rate(hist: dict) -> float:
    num_kmers_appear_once = hist.get(2)
    num_unique_kmers = find_unique_kmers(hist)
    error_rate = num_kmers_appear_once/num_unique_kmers

    if error_rate == 0:
        error_rate = 1

    return error_rate


def find_unique_kmers(hist: dict) -> int:
    return sum(hist.values())


def apply_savgol_filt(input: List[float]) -> List[float]:
    # Questionable Values, but experimental nevertheless.
    return savgol_filter(input, 3, 2)


def find_maxima(input_values: List[int]) -> int:
    # Gives you a list of indexs where contains a maximum.
    list_of_maximas = argrelextrema(input_values, np.greater)[0]

    # Place holder value for the current highest value.
    curr_largest_maxima = 0

    for curr_maxima in list_of_maximas:
        x = input_values[curr_maxima]
        if curr_largest_maxima < x:
            curr_largest_maxima = x

    return int(curr_largest_maxima)
