import numpy as np
from scipy.stats import binom
    

def compute_weight(k, num_hashes, p_val = 0.01, mut_thresh = 0.05, coverage = 1):
    """
    Computes the correct false-positive weight given various parameters.
    :param ksize: kmer size
    :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
    :param p_val: target quantile. If p_val <= 0, use expected unmutated kmers instead.
    :param mut_thresh: mutation probability threshold for species equivalence.
    :return: weight w, p_val quantile for number of non-mutated kmers.
    """
    non_mut_p = (1-mut_thresh)**k
    non_mut_quantile = binom.ppf(p_val, num_hashes, non_mut_p)
    non_mut_quantile_coverage = float(int(non_mut_quantile * coverage))
    w = non_mut_quantile_coverage / (num_hashes - non_mut_quantile_coverage)
    return w, non_mut_quantile, non_mut_quantile_coverage