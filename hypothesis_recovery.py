import numpy as np
import warnings
from scipy.stats import binom
from scipy.special import betaincinv
import pandas as pd
warnings.filterwarnings("ignore")


def get_nontrivial_idx(A, y):
    """
    This function takes the sparse matrix A and the vector y and returns the indices of the columns of A
    that have a non-zero inner product/overlap with y.
    :param A: sparse matrix
    :param y: vector
    :return: indices of columns of A that have a non-zero inner product with y
    """
    inners = A.T @ y
    nonz_idx = np.nonzero(inners)[0]
    return nonz_idx


def get_exclusive_indicators(A):
    """
    This function takes the sparse matrix A and returns a list of lists,
    where each of the ith list is the set of rows that are non-zero in the ith column.
    :param A: A sparse matrix. Should be binary, but doesn't have to be.
    :return: list(list(int))
    """
    # TODO: currently this isn't finding the unique rows, but rather the rows that are non-zero
    unique_locs = []
    m, N = A.shape
    # sum all the columns up
    col_sums = A.sum(axis=1)
    # look for the rows that have a 1 in them
    #unique_rows = np.nonzero(col_sums > 0)[0]  # FIXME: this is looking for rows that have an entry > 0, not == 1
    unique_rows = np.nonzero(col_sums == 1)[0]
    # turn this into a set
    unique_rows = set(unique_rows)
    # for each column, find the rows that are non-zero
    for i in range(N):
        non_zero_locs = np.nonzero(A[:, i])[0]
        # find the intersection of the two sets
        unique_in_col = list(unique_rows.intersection(non_zero_locs))
        # sort this, if need be
        # unique_in_col = sorted(unique_in_col)
        unique_locs.append(unique_in_col)
    return unique_locs


def get_alt_mut_rate(nu, thresh, ksize, significance=0.99):
    """
    Computes the alternative mutation rate for a given significance level. I.e. how much higher would the mutation rate
    have needed to be in order to have a false positive rate of significance (since we are setting the false negative
    rate to significance by design)?
    :param nu: Number of k-mers exclusive to the organism under consideration
    :param thresh: Number of exclusive k-mers I would need to observe in order to reject the null hypothesis (i.e.
    accept that the organism is present)
    :param ksize: K-mer size
    :param significance: value between 0 and 1 expressing the desired false positive rate (and by design, the false
    negative rate)
    :return: float (alternative mutation rate; how much higher would the mutation rate have needed to be in order to
    make FP and FN rates equal to significance)
    """
    # Replace binary search with the regularized incomplete Gamma function inverse: Solve[significance ==
    #   BetaRegularized[1 - (1 - mutCurr)^k, nu - thresh,
    #    1 + thresh], mutCurr]
    # per mathematica
    mut = 1 - (1 - betaincinv(nu - thresh, 1 + thresh, significance))**(1/ksize)
    if np.isnan(mut):
        return -1
    else:
        return mut


def single_hyp_test(
    A,
    y,
    unique_idx,
    ksize,
    significance=0.99,
    ani_thresh=0.95,
    min_coverage=1
):
    nu = len(unique_idx)
    
    non_mut_p = (ani_thresh)**ksize
    non_mut_thresh = binom.ppf(1-significance, nu, non_mut_p)
    act_conf = 1-binom.cdf(non_mut_thresh, nu, non_mut_p)
    nu_coverage = int(nu * min_coverage)
    non_mut_thresh_coverage = binom.ppf(1-significance, nu_coverage, non_mut_p)
    act_conf_coverage = 1-binom.cdf(non_mut_thresh_coverage, nu_coverage, non_mut_p)
    
    alt_mut = get_alt_mut_rate(nu, non_mut_thresh, ksize, significance=significance)
    alt_mut_cover = get_alt_mut_rate(nu_coverage, non_mut_thresh_coverage, ksize, significance=significance)
    
    num_matches = len(np.nonzero(y[unique_idx])[0])
    p_val = binom.cdf(num_matches, nu, non_mut_p)
    is_present = (num_matches >= non_mut_thresh_coverage)
    return is_present, p_val, nu, nu_coverage, num_matches, non_mut_thresh, non_mut_thresh_coverage, act_conf, act_conf_coverage, alt_mut, alt_mut_cover


def hypothesis_recovery(
    A,
    y,
    ksize,
    significance=0.99,
    ani_thresh=0.95,
    min_coverage=1,
):
    nont_idx = get_nontrivial_idx(A, y)
    N = np.shape(A)[1]
    A_sub = A[:,nont_idx]
    exclusive_indicators = get_exclusive_indicators(A_sub)
    
    nontriv_flags = np.zeros(N)
    nontriv_flags[nont_idx] = 1
    
    is_present = np.zeros(N)
    p_vals = np.zeros(N)
    alt_probs = np.zeros(N)
    num_unique_kmers = np.zeros(N)
    num_unique_kmers_coverage = np.zeros(N)
    num_matches = np.zeros(N)
    raw_thresholds = np.zeros(N)
    coverage_thresholds = np.zeros(N)
    act_conf = np.zeros(N)
    act_conf_coverage = np.zeros(N)
    alt_mut = np.zeros(N)
    alt_mut_cover = np.zeros(N)
    
    for i in range(len(nont_idx)):
        exclusive_idx = exclusive_indicators[i]
        curr_result = single_hyp_test(
            A_sub,
            y,
            exclusive_idx,
            ksize,
            significance=significance,
            ani_thresh=ani_thresh,
            min_coverage=min_coverage,
        )
        curr_idx = nont_idx[i]
        
        is_present[curr_idx], p_vals[curr_idx], num_unique_kmers[curr_idx], num_unique_kmers_coverage[curr_idx], num_matches[curr_idx], raw_thresholds[curr_idx], coverage_thresholds[curr_idx], act_conf[curr_idx], act_conf_coverage[curr_idx], alt_mut[curr_idx], alt_mut_cover[curr_idx] = curr_result
    
    return is_present, p_vals, num_unique_kmers, num_unique_kmers_coverage, num_matches, raw_thresholds, coverage_thresholds, act_conf, act_conf_coverage, alt_mut, alt_mut_cover, nontriv_flags


def hypothesis_recovery_pandas(
    A,
    y,
    ksize,
    significance=0.99,
    ani_thresh=0.95,
    min_coverage=1,
):
    # This will be the same as the above function, but using a pandas dataframe instead of using a bunch of numpy arrays and named variables
    nont_idx = get_nontrivial_idx(A, y)
    N = np.shape(A)[1]
    A_sub = A[:, nont_idx]
    exclusive_indicators = get_exclusive_indicators(A_sub)
    nontriv_flags = np.zeros(N)
    nontriv_flags[nont_idx] = 1
    # Create a pandas dataframe to store the results
    results = pd.DataFrame(
        index=nont_idx,
        columns=[
            'is_present',
            'p_vals',
            'num_unique_kmers',
            'num_unique_kmers_coverage',
            'num_matches',
            'raw_thresholds',
            'coverage_thresholds',
            'act_conf',
            'act_conf_coverage',
            'alt_mut',
            'alt_mut_cover',
        ]
    )
    for i in range(len(nont_idx)):
        exclusive_idx = exclusive_indicators[i]
        curr_result = single_hyp_test(
            A_sub,
            y,
            exclusive_idx,
            ksize,
            significance=significance,
            ani_thresh=ani_thresh,
            min_coverage=min_coverage,
        )
        curr_idx = nont_idx[i]
        results.loc[curr_idx] = curr_result
    return results, nontriv_flags