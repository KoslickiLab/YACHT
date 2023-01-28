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


# def recover_abundance_from_vectors(A, y, w):
#     """
#     Runs the linear program for quantile regression with weight w on the equation Ax = y.
#     :param A: matrix (reference database)
#     :param y: vector (sample kmer counts)
#     :param w: False positive weight
#     :return: vector x (estimated organism counts)
#     """
#     K, N = np.shape(A)
#     x = cp.Variable(N)
#     u = cp.Variable(K)
#     v = cp.Variable(K)
#     tau = 1 / (w + 1)
#     ones_K = np.ones(K)
#     objective = cp.Minimize(
#         tau * (ones_K @ u) + (1 - tau) * (ones_K @ v)
#     )
#     constraints = [
#         x >= 0,
#         u >= 0,
#         v >= 0,
#         u - v + (A @ x) == y,
#     ]
#     prob = cp.Problem(objective, constraints)
#     result = prob.solve(solver=cp.SCIPY, verbose=False)
#     recov_y = A @ x.value
#     resid = y - (A @ x.value)
#     return x.value, resid
    

# def load_reference_metadata(
#     matrix_file,
#     ksize,
# ):
#     prefix = matrix_file.split('ref_matrix_processed.npz')[0]
#     hash_to_idx_file = prefix + 'hash_to_col_idx.csv'
#     processed_org_file = prefix + 'processed_org_idx.csv'
    
#     reference_matrix = load_npz(matrix_file)
#     hash_to_idx = utils.load_hashes(hash_to_idx_file)
#     organism_data = pd.read_csv(processed_org_file)
    
#     return reference_matrix, hash_to_idx, hash_to_idx_file, organism_data


# def recover_abundance_data(
#     ref_matrix,
#     sample_vector,
#     ref_organism_data,
#     ksize,
#     mut_thresh,
#     significance,
#     num_kmers_quantile,
#     min_coverage,
#     num_sample_kmers,
#     num_unique_sample_kmers,
#     sample_scale,
#     w=None,
# ):
#     recov_org_data = ref_organism_data.copy()
#     recov_org_data['num_total_kmers_in_sample_sketch'] = num_sample_kmers
#     recov_org_data['num_unique_kmers_in_sample_sketch'] = num_unique_sample_kmers
#     recov_org_data['sample_scale_factor'] = sample_scale
#     #recov_org_data['num_total_kmers_in_sample_sketch_scaled'] = num_sample_kmers*sample_scale
    
#     sample_diff_idx = np.nonzero(np.array(np.abs(recov_org_data['sample_scale_factor'] - recov_org_data['genome_scale_factor'])))[0]
#     sample_diffs = list(recov_org_data['organism_name'][sample_diff_idx])
    
#     if len(sample_diffs) > 0:
#         raise ValueError('Sample scale factor does not equal genome scale factor for organism %s and %d others.'%(sample_diffs[0],len(sample_diffs)-1))
    
#     est_count_genomes = np.round(num_sample_kmers / np.mean(recov_org_data['num_total_kmers_in_genome_sketch']))
#     recov_org_data['est_count_genomes_in_sample'] = est_count_genomes
    
#     if w is None:
#         num_kmers_for_pval = int(np.quantile(recov_org_data['num_unique_kmers_in_genome_sketch'], num_kmers_quantile))
#         recov_org_data['num_unique_kmers_for_pval'] = num_kmers_for_pval
#         w, min_quantile, min_quantile_cov = cw.compute_weight(ksize, num_kmers_for_pval, p_val = p_val, mut_thresh = mut_thresh, coverage = min_coverage)
#         recov_org_data['unmutated_kmer_quantile'] = min_quantile
#         recov_org_data['unmutated_kmer_quantile_with_coverage'] = min_quantile_cov
#     else:
#         warnings.warn('w set manually; specified p_val overriden.')
#         recov_org_data['num_unique_kmers_for_pval'] = -1
#         recov_org_data['unmutated_kmer_quantile'] = -1
#         recov_org_data['unmutated_kmer_quantile_with_coverage'] = -1
    
#     recov_org_data['min_coverage'] = min_coverage
#     recov_org_data['w'] = w
    
#     abundance, residual = recover_abundance_from_vectors(ref_matrix, sample_vector, w)  
#     recov_org_data['recovered_kmer_abundance'] = abundance    
#     recov_org_data['recovered_count_abundance'] = abundance/recov_org_data['num_total_kmers_in_genome_sketch']
    
#     recov_sample = ref_matrix @ recov_org_data['recovered_kmer_abundance']
    
#     sample_nonzero = np.nonzero(sample_vector)[0]
# #     #overestimates correspond to mutations
#     overestimates = np.maximum(recov_sample - sample_vector, 0)
# #     #underestimates correspond to missed kmers
#     underestimates = np.maximum(sample_vector - recov_sample, 0)
# #     #we count underestimates where kmers are missed entirely:
#     under_non_recov = underestimates[recov_sample == 0]
    
#     recov_org_data['total_sample_kmers_in_ref'] = np.sum(sample_vector)
#     recov_org_data['recovery_sample_overestimates'] = np.sum(overestimates)
#     recov_org_data['recovery_sample_overestimates'] = np.sum(underestimates)
#     recov_org_data['recovery_sample_missed_kmers'] = np.sum(under_non_recov)
    
#     recov_org_data['est_mut_kmers_in_sample'] = recov_org_data['recovery_sample_overestimates']/recov_org_data['sample_scale_factor']
#     recov_org_data['est_known_kmers_in_sample'] =  recov_org_data['total_sample_kmers_in_ref'] - recov_org_data['recovery_sample_missed_kmers'] + recov_org_data['est_mut_kmers_in_sample']
    
#     recov_org_data['recovery_unknown_pct_est'] = 1 - recov_org_data['est_known_kmers_in_sample']/    recov_org_data['num_total_kmers_in_sample_sketch']
    
#     return recov_org_data, abundance, recov_sample, overestimates, underestimates


# def recover_abundance_from_files(
#     matrix_file,
#     sample_file,
#     ksize,
#     mut_thresh,
#     significance,
#     num_kmers_quantile,
#     min_coverage,
#     output_filename=None,
#     w=None
# ):
#     """
#     Runs linear program for unknown estimation off of files generated by ref_matrix.py and creates human-readable results file.
#     :param matrix_file: location of ref_matrix_processed.npz file (A matrix)
#     :param sample_file: location of sample.sig file (y vector)
#     :param ksize: kmer size
#     :param mut_thresh: mutation cutoff for species equivalence
#     :param significance: minimum probability of individual true negative.
#     :param num_kmers_quantile: quantile for determining representative number of kmers in sketch to be used in calculation of p-value.
#     :param output_filename: destination for results file; if blank, no file will be written
#     :param w: false positive weight. Optional; if set, overrides significance.
#     :return: pandas dataframe containing recovered abundances and metadata.
#     """
#     (
#         reference_matrix,
#         hash_to_idx,
#         hash_to_idx_file,
#         organism_data
#     ) = load_reference_metadata(matrix_file, ksize)
    
#     sample_vector, sample_sig, num_kmers_non_ref_unique, num_kmers_non_ref_total = sv.sample_vector_from_files(sample_file, hash_to_idx_file, ksize)
#     sample_scale = sample_sig.minhash.scaled
#     num_sample_kmers = utils.get_num_kmers(sample_sig, scale = False)
#     num_unique_sample_kmers = len(list(sample_sig.minhash.hashes))

#     recov_org_data, abundance, recov, over, under = recover_abundance_data(
#         reference_matrix,
#         sample_vector,
#         organism_data,
#         ksize,
#         mut_thresh,
#         significance,
#         num_kmers_quantile,
#         min_coverage,
#         num_sample_kmers,
#         num_unique_sample_kmers,
#         sample_scale,
#         w=w,
#     )
    
#     if output_filename:
#         recov_org_data.to_csv(output_filename)

#     return recov_org_data


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description="This script estimates the abundance of microorganisms from a reference database matrix and metagenomic sample.",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument('--ref_file', help='Reference database matrix in npz format', required=True)
#     parser.add_argument('--ksize', type=int, help='Size of kmers used in sketch', required=True)
#     parser.add_argument('--sample_file', help='Metagenomic sample in .sig format', required=True)
#     parser.add_argument('--w', type=float, help='False positive weight. If set manually, overrides significance argument.', required=False, default = None)
#     parser.add_argument('--mut_thresh', type=float, help='mutation cutoff for species equivalence.', required=False, default = 0.05)
#     parser.add_argument('--significance', type=float, help='Minimum probability of individual true negative.', required=False, default = 0.01)
#     parser.add_argument('--num_kmers_quantile', type=float, help='To compute false negative p-val, assume each organism has constant number of kmers in the sketch given by this quantile of the actual kmer counts.', required=False, default = 0.33)
#     parser.add_argument('--min_coverage', type=float, help='To compute false negative weight, assume each organism has this minimum coverage in sample. Should be between 0 and 1.', required=False, default = 1)
#     parser.add_argument('--outfile', help='csv destination for results', required=True)
#     args = parser.parse_args()
#     recover_abundance_from_files(
#         args.ref_file,
#         args.sample_file,
#         args.ksize,
#         args.mut_thresh,
#         args.significance,
#         args.num_kmers_quantile,
#         args.min_coverage,
#         args.outfile,
#         w = args.w
#     )
