from bio_hansel.kmer_val_calc.kmer_utils import calc_error_rate, calc_avg_kmer_depth
from scipy.stats import norm as stats


def find_min_kmer_val(hist: dict) -> float:
    error_rate = calc_error_rate(hist)
    kmer_cov_freq = calc_avg_kmer_depth(hist)
    min_kmer_val = stats.ppf(0.95, kmer_cov_freq * error_rate)
    return min_kmer_val
