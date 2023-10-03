def recover_abundance_from_vectors(A, y, w):
    """
    Runs the linear program for quantile regression with weight w on the equation Ax = y.
    :param A: matrix (reference database)
    :param y: vector (sample kmer counts)
    :param w: False positive weight
    :return: vector x (estimated organism counts)
    """
    K, N = np.shape(A)
    x = cp.Variable(N)
    u = cp.Variable(K)
    v = cp.Variable(K)
    tau = 1 / (w + 1)
    ones_K = np.ones(K)
    objective = cp.Minimize(
        tau * (ones_K @ u) + (1 - tau) * (ones_K @ v)
    )
    constraints = [
        x >= 0,
        u >= 0,
        v >= 0,
        u - v + (A @ x) == y,
    ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.SCIPY, verbose=False)
    recov_y = A @ x.value
    resid = y - (A @ x.value)
    return x.value, resid


import numpy as np
import sourmash
import utils
import pdb


def sample_vector_from_signature(signature, hash_to_idx, normalize=False):
    """
    Given a signature and a dictionary mapping hashes to indices, return a vector whose basis is indexed by
    all hashes seen in the training dictionary. The value at each index is the abundance of the hash in the
    signature.

    :param signature: a sourmash minhash object
    :param hash_to_idx: a dictionary mapping hashes of the training dictionary to indices
    :return: numpy vector
    """
    K = len(list(hash_to_idx.keys()))
    sample_vec = np.zeros(K)
    sig_hashes = signature.minhash.hashes
    sig_hash_overlap = np.intersect1d(sig_hashes, list(hash_to_idx.keys()))
    sig_hash_diff = np.setdiff1d(sig_hashes, sig_hash_overlap)
    sig_hash_diff_values = [sig_hashes[h] for h in sig_hash_diff]
    num_hash_diff_unique = len(list(sig_hash_diff))
    num_hash_diff_total = np.sum(sig_hash_diff_values)
    for sh in sig_hash_overlap:
        idx = hash_to_idx[sh]
        sample_vec[idx] = sig_hashes[sh]
    if normalize:
        sample_vec = sample_vec / utils.get_num_kmers(signature, scale = False)
    return sample_vec, num_hash_diff_unique, num_hash_diff_total


def sample_vector_from_files(sig_filename, hash_filename, ksize):
    """
    Helper function to load a signature and a hash_to_idx dictionary from files and return a sample vector.

    :param sig_filename: filename of the sourmash signature
    :param hash_filename: filename of the hash_to_col_idx.csv file which maps hashes to indices
    :param ksize: ksize of the signature
    :return: numpy vector (sample vector y)
    """
    sample_sig = utils.load_signature_with_ksize(sig_filename, ksize)
    hash_to_idx = utils.load_hashes_to_index(hash_filename)
    sample_vector, num_hash_diff_unique, num_hash_diff_total = sample_vector_from_signature(sample_sig, hash_to_idx)
    return sample_vector, sample_sig, num_hash_diff_unique, num_hash_diff_total


# This was in compute_weight.py
import numpy as np
from scipy.stats import binom


def compute_weight(k, num_hashes, p_val=0.01, mut_thresh=0.05, coverage=1):
    """
    Computes the correct false-positive weight given various parameters.
    :param ksize: kmer size
    :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
    :param p_val: target quantile. If p_val <= 0, use expected unmutated kmers instead.
    :param mut_thresh: mutation probability threshold for species equivalence.
    :return: weight w, p_val quantile for number of non-mutated kmers.
    """
    non_mut_p = (1 - mut_thresh) ** k
    non_mut_quantile = binom.ppf(p_val, num_hashes, non_mut_p)
    non_mut_quantile_coverage = float(int(non_mut_quantile * coverage))
    w = non_mut_quantile_coverage / (num_hashes - non_mut_quantile_coverage)
    return w, non_mut_quantile, non_mut_quantile_coverage

# This was in solve_lp.py
import numpy as np
import cvxpy as cp
import pandas as pd
import csv
import compute_weight as cw
from scipy.sparse import load_npz
import argparse
import warnings


def solve_lp(A, y, w):
    """
    Runs the linear program for nonnegative quantile regression with weight w on the equation Ax = y.
    :param A: matrix (reference database)
    :param y: vector (sample kmer counts)
    :param w: False positive weight vector
    :return: vector x (estimated organism kmer counts)
    """
    K, N = np.shape(A)
    x = cp.Variable(N)
    u = cp.Variable(K)
    v = cp.Variable(K)
    tau = 1 / (w + 1)
    objective = cp.Minimize(
        tau @ u + (1 - tau) @ v
    )
    constraints = [
        x >= 0,
        u >= 0,
        v >= 0,
        u - v + (A @ x) == y,
    ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.SCIPY, verbose=False)
    recov_y = A @ x.value
    resid = y - (A @ x.value)
    return x.value, resid