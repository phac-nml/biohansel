from bio_hansel.kmer_val_calc.kmer_utils import calc_error_rate, calc_avg_kmer_depth
from scipy.stats import norm as stats


def find_min_kmer_val(hist: dict) -> float:
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

    :returns:
            The minimum kmer cut off value.
    """
    error_rate = calc_error_rate(hist)
    kmer_cov_freq = calc_avg_kmer_depth(hist)
    min_kmer_val = stats.ppf(0.95, kmer_cov_freq * error_rate)
    return min_kmer_val
