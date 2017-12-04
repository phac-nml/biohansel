from scipy.signal import argrelextrema, savgol_filter
from typing import List
import numpy as np

from bio_hansel.subtyping_params import SubtypingParams


def calc_avg_kmer_depth(hist: dict, subtyping_params: SubtypingParams) -> int:
    """Calculate the average k-mer depth
    Note:
            To find the average kmer coverage value, we use the assumption that most data follows the following curve:
                |.
                | .       |------- Here is the average k-mer coverage depth value.
                | .     . | .
                | .    .  |  .
                | .   .   |   .
                |  . .    |    . .
                |   .     |       . . .
                |______________________________

            Where the maxima's value would be your average k-mer coverage value by accessing it's key.

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.
            :param subtyping_params: Object containing parameters for savgol filter

    :returns:
            The average k-mer coverage depth.
    """
    maximum = find_maxima(apply_savgol_filt(list(hist.values()), subtyping_params))
    kmer_coverage = -1

    for key, value in hist.items():
        if maximum == value:
            kmer_coverage = key
            break

    return kmer_coverage


def calc_error_rate(hist: dict) -> float:
    """Calculates the estimated error rate through observation.
    Note:
            As we assume that k-mers that appear once are errors, we simply calculate the error rate by
            kmers that appear once / unique kmers.

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.

    :return:
            float value containing the error rate.
    """
    num_kmers_appear_once = hist.get(2)
    num_unique_kmers = find_unique_kmers(hist)
    error_rate = num_kmers_appear_once/num_unique_kmers

    # round up if there's a 0 error rate.
    if error_rate == 0:
        error_rate = 1

    return error_rate


def find_unique_kmers(hist: dict) -> int:
    """Find the number of unique k-mers within the observation
    Note:
            We can simply follow the "curve of the graph" to find the number of unique k-mers.

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.

    :return:
            int containing the number of unique kmers
    """
    return sum(hist.values())


def apply_savgol_filt(input: List[float], subtyping_params: SubtypingParams) -> np.ndarray:
    """Apply savgol filter to the observation to smooth out curve.

    :arg:
            :param input: List containing the values of the y-axis of the k-mer frequency graph.
            :param subtyping_params: Parameters for the savgol filter method.

    :return:
            List containing the smoothed values of the y-axis of the k-mer frequency graph.
    """
    # Savgol filter takes the input values, window length, and degree of polynomial.
    return savgol_filter(input, subtyping_params.savgol_window_len, subtyping_params.savgol_poly_degree)


def find_maxima(input_values: np.ndarray) -> int:
    """Finds the local maxima containing the average k-mer coverage depth value.

    :arg:
            :param input_values: List containing the values of the y-axis of the k-mer frequency graph.

    :return:
            average k-mer coverage depth value
    """
    # Gives you a list of indexs where contains a maximum.
    list_of_maximas = argrelextrema(input_values, np.greater)[0]

    # Place holder value for the current highest value.
    curr_largest_maxima = 0

    for curr_maxima in list_of_maximas:
        x = input_values[curr_maxima]
        if curr_largest_maxima < x:
            curr_largest_maxima = x

    return int(curr_largest_maxima)
