import logging

import numpy as np
import pandas as pd
from scipy.stats import poisson as poisson
from scipy.signal import argrelextrema, savgol_filter

from ..subtyping_params import SubtypingParams


def mean_kmer_depth(hist: pd.DataFrame, subtyping_params: SubtypingParams) -> int:
    """Calculate the mean k-mer coverage depth

    Note:
            To find the mean kmer coverage value, we use the assumption that most data follows the following curve:

                |.
                | .       |------- Here is the mean k-mer coverage depth value.
                | .     . | .
                | .    .  |  .
                | .   .   |   .
                |  . .    |    . .
                |   .     |       . . .
                |______________________________

            The maximum value is the mean k-mer coverage value by accessing it's key.

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.
            :param subtyping_params: Object containing parameters for savgol filter

    :returns:
            The average k-mer coverage depth.
    """
    # Apply savgol filter to the observation to smooth out curve.
    savgol_out = savgol_filter(hist.obs,
                               window_length=subtyping_params.savgol_window_len,
                               polyorder=subtyping_params.savgol_poly_degree)
    logging.debug('savgol_out %s', savgol_out)
    maximum = find_mean_kmer_cov(savgol_out)
    logging.debug('maximum in savgol %s', maximum)
    return int(hist[hist.obs == maximum].freq)


def error_rate(hist: pd.DataFrame) -> float:
    """Calculates the estimated k-mer error rate through observation.
    Note:
            As we assume that k-mers that appear once are errors, we simply calculate the error rate by
            kmers that appear once / unique kmers.

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.

    :return:
            float value containing the error rate.
    """
    n_singleton_kmers = float(hist[hist.freq == 1].obs)
    logging.debug('n_singleton_kmers %s', n_singleton_kmers)
    # get the number of unique k-mers within the observation
    #   We can simply follow the "curve of the graph" to find the number of unique k-mers.
    n_unique_kmers = (hist['freq'] * hist['obs']).sum()
    logging.debug('n_unique_kmers %s', n_unique_kmers)
    # case by / 0
    try:
        e = n_singleton_kmers / n_unique_kmers
    except ZeroDivisionError:
        logging.warning("Could not calculate error rate. {} = n_singleton_kmers, {} = n_unique_kmers"
                        .format(n_singleton_kmers, n_unique_kmers))
        raise ValueError('Could not calculate error_rate.')

    logging.debug('error_rate %s', e)
    # round up if there's a 0 error rate.

    return e


def find_mean_kmer_cov(xs: np.ndarray) -> int:
    """Finds the local maxima containing the average k-mer coverage depth value.

    :arg:
            :param xs: List containing the values of the y-axis of the k-mer frequency graph.

    :return:
            average k-mer coverage depth value
    """
    # Gives you a list of indexs where contains a maximum.

    i_xs = argrelextrema(xs, np.greater)
    logging.debug('argrelextrema %s', i_xs)
    if isinstance(i_xs, tuple):
        i_xs = i_xs[0]
    m = int(np.max(xs[i_xs]))
    logging.debug('max %s', m)
    return m


def min_kmer_freq_threshold(hist: pd.DataFrame, subtyping_params: SubtypingParams) -> float:
    """Find the minimum kmer frequency coverage value from observation
    Note:
            This method will find three things:
            1) Error Rate
            2) The average k-mer coverage depth
            3) The minimum kmer cut off value

            How does this work?
            1) Estimate the expected k-mer coverage depth from observation.
            2) (Slightly under-) estimate the k-mer error rate by assuming k-mers that appear once are errors.
            3) Multiply those two things together to get the (slightly underestimated)
            expected k-mer coverage depth of errors.
            4) Plug the expected k-mer coverage depth of errors into Poisson to pull out the minimum depth
            to be confident "the observation is not entirely cause by errors".

    :arg:
            :param hist: Dictionary containing the histogram of the observed kmer values.
            :param subtyping_params: Parameters containing savgol and percent confidence values.

    :returns:
            The minimum kmer cut off value.
    """
    e = error_rate(hist)
    logging.debug('error rate: %s', e)
    d = mean_kmer_depth(hist, subtyping_params)
    logging.debug('kmer_cov_freq: %s', d)
    min_kmer_val = poisson.ppf(subtyping_params.kmer_cov_perc_confidence, d * e)
    logging.debug('min_kmer_val %s', min_kmer_val)
    return min_kmer_val