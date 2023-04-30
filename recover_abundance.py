import numpy as np
import pandas as pd
import hypothesis_recovery as hr
from scipy.sparse import load_npz
import argparse
import utils
import warnings
import os
warnings.filterwarnings("ignore")


def recover_abundance_data_hyp(
    ref_matrix,
    sample_vector,
    ref_organism_data,
    ksize,
    ani_thresh,
    significance,
    min_coverage,
    num_sample_kmers,
    num_unique_sample_kmers,
    sample_scale,
):
    recov_org_data = ref_organism_data.copy()
    recov_org_data['num_total_kmers_in_sample_sketch'] = num_sample_kmers
    recov_org_data['num_unique_kmers_in_sample_sketch'] = num_unique_sample_kmers
    recov_org_data['sample_scale_factor'] = sample_scale
    
    sample_diff_idx = np.nonzero(np.array(np.abs(recov_org_data['sample_scale_factor'] - recov_org_data['genome_scale_factor'])))[0]
    sample_diffs = list(recov_org_data['organism_name'][sample_diff_idx])
    if len(sample_diffs) > 0:
        raise ValueError('Sample scale factor does not equal genome scale factor for organism %s and %d others.'%(sample_diffs[0],len(sample_diffs)-1))
    
    recov_org_data['min_coverage'] = min_coverage
    
    in_sample_est, p_vals, nu, nu_coverage, num_matches, acceptance_threshold_wo_coverage, acceptance_threshold_with_coverage, actual_confidence_wo_coverage, actual_confidence_with_coverage, alt_mut, alt_confidence_mut_rate_with_coverage, nontriv_flags = hr.hypothesis_recovery(ref_matrix, sample_vector, ksize, significance=significance, ani_thresh=ani_thresh, min_coverage=min_coverage)
    
    recov_org_data['nontrivial_overlap'] = nontriv_flags
    recov_org_data['in_sample_est'] = in_sample_est
    recov_org_data['num_exclusive_kmers'] = nu
    recov_org_data['num_exclusive_kmers_with_coverage'] = nu_coverage
    recov_org_data['num_matches'] = num_matches
    recov_org_data['acceptance_threshold_wo_coverage'] = acceptance_threshold_wo_coverage
    recov_org_data['acceptance_threshold_with_coverage'] = acceptance_threshold_with_coverage
    recov_org_data['actual_confidence_wo_coverage'] = actual_confidence_wo_coverage
    recov_org_data['actual_confidence_with_coverage'] = actual_confidence_with_coverage
    recov_org_data['p_vals'] = p_vals
    recov_org_data['alt_confidence_mut_rate'] = alt_mut
    recov_org_data['alt_confidence_mut_rate_with_coverage'] = alt_confidence_mut_rate_with_coverage
    return recov_org_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script estimates the abundance of microorganisms from a reference database matrix and metagenomic sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_matrix', help='Reference database matrix in npz format', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers used in sketch', required=True)
    parser.add_argument('--sample_file', help='Metagenomic sample in .sig format', required=True)
    parser.add_argument('--ani_thresh', type=float, help='mutation cutoff for species equivalence.',
                        required=False, default=0.95)
    parser.add_argument('--significance', type=float, help='Minimum probability of individual true negative.',
                        required=False, default=0.99)
    parser.add_argument('--min_coverage', type=float, help='To compute false negative weight, assume each organism has this minimum coverage in sample. Should be between 0 and 1.', required=False, default = 1)
    parser.add_argument('--outfile', help='csv destination for results', required=True)

    # parse the arguments
    args = parser.parse_args()
    ref_matrix = args.ref_matrix  # location of ref_matrix_processed.npz file (A matrix)
    sample_file = args.sample_file  # location of sample.sig file (y vector)
    ksize = args.ksize
    ani_thresh = args.ani_thresh  # ANI cutoff for species equivalence
    significance = args.significance  # Minimum probability of individual true negative.
    min_coverage = args.min_coverage  # assume each organism has this minimum coverage in sample. Should be between 0 and 1
    outfile = args.outfile  # csv destination for results

    # check that ksize is an integer
    if not isinstance(ksize, int):
        raise ValueError('ksize must be an integer.')
    # check if min_coverage is between 0 and 1
    if min_coverage < 0 or min_coverage > 1:
        raise ValueError('min_coverage must be between 0 and 1.')

    # Get the training data names
    prefix = ref_matrix.split('ref_matrix_processed.npz')[0]
    hash_to_idx_file = prefix + 'hash_to_col_idx.csv'
    processed_org_file = prefix + 'processed_org_idx.csv'

    # make sure all these files exist
    if not os.path.exists(ref_matrix):
        raise ValueError(f'Reference matrix file {ref_matrix} does not exist. Please run ref_matrix.py first.')
    if not os.path.exists(hash_to_idx_file):
        raise ValueError(f'Hash to index file {hash_to_idx_file} does not exist. Please run ref_matrix.py first.')
    if not os.path.exists(processed_org_file):
        raise ValueError(
            f'Processed organism file {processed_org_file} does not exist. Please run ref_matrix.py first.')

    # load the training data
    reference_matrix = load_npz(ref_matrix)
    hash_to_idx = utils.load_hashes(hash_to_idx_file)
    organism_data = pd.read_csv(processed_org_file)

    # get the sample y vector (indexed by hash/k-mer, with entry = number of times k-mer appears in sample)
    sample_sig = utils.load_signature_with_ksize(sample_file, ksize)
    # total number of hashes in the training dictionary
    K = len(list(hash_to_idx.keys()))
    # initialize the sample vector
    sample_vector = np.zeros(K)
    # get the hashes in the signature (it's for a single sample)
    sample_hashes = sample_sig.minhash.hashes
    # get the hashes that are in both the sample and the training dictionary
    sample_intersect_training_hashes = np.intersect1d(sample_hashes, list(hash_to_idx.keys()))
    for sh in sample_intersect_training_hashes:
        idx = hash_to_idx[sh]
        sample_vector[idx] = sample_hashes[sh]

    # get the number of kmers in the sample from the scaled sketch
    sample_scale = sample_sig.minhash.scaled
    num_sample_kmers = utils.get_num_kmers(sample_sig, scale=False)
    # get the number of unique kmers in the sample
    num_unique_sample_kmers = len(list(sample_sig.minhash.hashes))

    # prep the output data structure, copying over the organism data
    recov_org_data = organism_data.copy()
    recov_org_data['num_total_kmers_in_sample_sketch'] = num_sample_kmers
    recov_org_data['num_unique_kmers_in_sample_sketch'] = num_unique_sample_kmers
    recov_org_data['sample_scale_factor'] = sample_scale

    # check that the sample scale factor is the same as the genome scale factor for all organisms
    sample_diff_idx = \
        np.nonzero(np.array(np.abs(recov_org_data['sample_scale_factor'] - recov_org_data['genome_scale_factor'])))[0]
    sample_diffs = list(recov_org_data['organism_name'][sample_diff_idx])
    if len(sample_diffs) > 0:
        raise ValueError('Sample scale factor does not equal genome scale factor for organism %s and %d others.' % (
        sample_diffs[0], len(sample_diffs) - 1))

    recov_org_data['min_coverage'] = min_coverage

    hyp_recovery_df, nontriv_flags = hr.hypothesis_recovery(
        reference_matrix, sample_vector, ksize, significance=significance, ani_thresh=ani_thresh, min_coverage=min_coverage)

    # Boolean indicating whether genome shares at least one k-mer with sample
    recov_org_data['nontrivial_overlap'] = nontriv_flags

    # Main output: Boolean indicating whether genome is present in sample
    recov_org_data['in_sample_est'] = hyp_recovery_df['in_sample_est']

    # Number of k-mers exclusive to genome
    recov_org_data['num_exclusive_kmers'] = hyp_recovery_df['num_unique_kmers']

    # Number of k-mers exclusive to genome, multiplied by min_coverage parameter
    recov_org_data['num_exclusive_kmers_with_coverage'] = hyp_recovery_df['num_unique_kmers_coverage']

    # Size of intersection between exclusive k-mers and sample
    recov_org_data['num_matches'] = hyp_recovery_df['num_matches']

    # Acceptance threshold without adjusting for coverage
    recov_org_data['acceptance_threshold_wo_coverage'] = hyp_recovery_df['acceptance_threshold_wo_coverage']

    # Acceptance threshold with adjusting for coverage
    recov_org_data['acceptance_threshold_with_coverage'] = hyp_recovery_df['acceptance_threshold_with_coverage']

    # computed confidence without adjusting for coverage
    recov_org_data['actual_confidence_wo_coverage'] = hyp_recovery_df['actual_confidence_wo_coverage']

    # computed confidence with adjusting for coverage
    recov_org_data['actual_confidence_with_coverage'] = hyp_recovery_df['actual_confidence_with_coverage']

    # Probability of observing this or more extreme result at ANI threshold.
    recov_org_data['p_vals'] = hyp_recovery_df['p_vals']

    # Mutation rate such that at this mutation rate, false positive rate = p_val. Does not account for
    # min_coverage parameter.
    recov_org_data['alt_confidence_mut_rate'] = hyp_recovery_df['alt_confidence_mut_rate']

    # Mutation rate such that at this mutation rate, false positive rate = p_val, accounting for min_coverage parameter.
    recov_org_data['alt_confidence_mut_rate_with_coverage'] = hyp_recovery_df['alt_confidence_mut_rate_with_coverage']

    # save the results
    recov_org_data.to_csv(outfile)
