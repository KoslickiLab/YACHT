import numpy as np
from scipy.stats import binom

#crude monte carlo method
#keeping around for now for testing purposes
def estimate_w_test(k, num_hashes, mut_thresh=0.05, est_n_orgs=1000, p_val=0.01, n_tests=10000):
    prob = (1 - mut_thresh) ** k
    b = []
    for i in range(n_tests):
        b.append(min(np.random.binomial(num_hashes, prob, (est_n_orgs, 1))))
    min_est = np.quantile(b, p_val)
    w = min_est / (num_hashes - min_est)
    return w, min_est


def deletion_cdf(x, ksize, num_hashes, mut_thresh = 0.05, est_num_genomes = 1000):
    """
    Computes the cdf for the maximum of est_num_genomes-many binomials with success probability 1-(1-mut_thresh)**k.
    :param x: input for cdf
    :param ksize: kmer size
    :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
    :param mut_thresh: mutation probability threshold for species equivalence.
    :param est_num_genomes: estimated number of different genomes in the sample.
    :return: probability the above random variable is bounded by x
    """
    mut_prob = 1 - (1 - mut_thresh) ** ksize
    single_deletion_cdf = binom.cdf(x, num_hashes, mut_prob)
    return single_deletion_cdf ** est_num_genomes


def find_mut_quantile(ksize, num_hashes, p_val = 0.01, mut_thresh = 0.05, est_num_genomes = 1000, tol = 1e-10):
    """
    Using binary search, finds the p_val-quantile for the random variable with CDF given by deletion_cdf. Since the distribution is not continuous, exact solutions are unlikely, so function returns the smallest quantile for q greater than or equal to p_val.
    :param ksize: kmer size
    :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
    :param p_val: target quantile
    :param mut_thresh: mutation probability threshold for species equivalence.
    :param est_num_genomes: estimated number of different genomes in the sample.
    :param tol: search terminates if sufficiently close solution is found.
    :return: approximate p_val-quantile of the random variable.
    """
    p_val_comp = 1 - p_val
    lower = 1
    upper = num_hashes
    x_curr= np.ceil(num_hashes/2)
    p_curr = -1
    count = 0
    while (np.abs(p_val_comp - p_curr) > tol):
        count += 1
        p_curr = deletion_cdf(x_curr, ksize, num_hashes, mut_thresh, est_num_genomes)
        if p_curr > p_val_comp:
            upper = x_curr
            x_curr = lower + np.ceil((upper - lower)/2)
        else:
            lower = x_curr
            x_curr = lower + np.ceil((upper - lower)/2)
        if upper - lower <= 1:
            x_curr = upper
            break
    return x_curr


def non_mut_mean(k, num_hashes, mut_thresh = 0.05):
    return num_hashes * (1-mut_thresh)**k


def compute_weight(k, num_hashes, p_val = 0.01, mut_thresh = 0.05, est_num_genomes = 1000, tol = 1e-10):
    """
    Computes the correct false-positive weight given various parameters.
    :param ksize: kmer size
    :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
    :param p_val: target quantile. If p_val <= 0, use expected unmutated kmers instead.
    :param mut_thresh: mutation probability threshold for species equivalence.
    :param est_num_genomes: estimated number of different genomes in the sample.
    :param tol: search terminates if sufficiently close solution is found.
    :return: weight w, 1-p_val quantile for number of non-mutated kmers.
    """
    if p_val > 0:
        mut_quantile = find_mut_quantile(k, num_hashes, p_val = p_val, mut_thresh = mut_thresh, est_num_genomes = est_num_genomes, tol = tol)
        non_mut_quantile = num_hashes - mut_quantile
    else:
        non_mut_quantile = non_mut_mean(k, num_hashes, mut_thresh = mut_thresh)
    w = non_mut_quantile / (num_hashes - non_mut_quantile)
    return w, non_mut_quantile


def unmut_quantiles(p_val, kmer_counts, k, mut_thresh=0.05):
    non_mut_prob = (1-mut_thresh)**k
    unmut_qs = np.zeros(np.shape(kmer_counts)[0])
    for (i, kmer_ct) in enumerate(kmer_counts):
        unmut_qs[i]=binom.ppf(p_val, kmer_ct, non_mut_prob)
    return unmut_qs

def get_weights(p_val, kmer_counts, k, mut_thresh=0.05):
    quantiles = unmut_quantiles(p_val, kmer_counts, k, mut_thresh = mut_thresh)
    return quantiles/(kmer_counts - quantiles)