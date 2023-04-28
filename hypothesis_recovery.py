import numpy as np
import cvxpy as cp
import pandas as pd
import csv
import sample_vector as sv
import compute_weight as cw
from scipy.sparse import load_npz
import argparse
import utils
import warnings
from scipy.stats import binom
import pdb
warnings.filterwarnings("ignore")


def get_nontrivial_idx(A, y):
    inners = A.T @ y
    nonz_idx = np.nonzero(inners)[0]
    return nonz_idx


def get_exclusive_indicators(A):
    """
    This function takes the sparse matrix A and returns a list of lists,
    where each the ith list is the set of rows that are non-zero in the ith column.
    :param A: A sparse matrix. Should be binary, but doesn't have to be.
    :return: list(list(int))
    """
    unique_locs = []
    m, N = A.shape
    # sum all the columns up
    col_sums = A.sum(axis=1)
    # look for the rows that have a 1 in them
    unique_rows = np.nonzero(col_sums > 0)[0]
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


def get_alt_mut_rate(nu, thresh, ksize, significance = 0.99, max_iters = 1000, epsi = 1e-10):
    upper = 1
    lower = 0
    prob = 1
    iters = 0
    while(np.abs(prob - significance) > epsi):
        mut_curr = (upper+lower)/2
        p_curr = (1-mut_curr)**ksize
        prob = binom.cdf(thresh, nu, p_curr)
        if prob > significance:
            upper = mut_curr
        else:
            lower = mut_curr
        iters += 1
        if iters > max_iters:
            return -1
    return mut_curr


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
