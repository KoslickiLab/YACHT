import numpy as np
import pandas as pd
import hypothesis_recovery as hr
from scipy.sparse import load_npz
import argparse
import utils
import warnings
import os
warnings.filterwarnings("ignore")


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
    hash_to_idx = utils.load_hashes(hash_filename)
    sample_vector, num_hash_diff_unique, num_hash_diff_total = sample_vector_from_signature(sample_sig, hash_to_idx)
    return sample_vector, sample_sig, num_hash_diff_unique, num_hash_diff_total


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
    
    is_present, p_vals, nu, nu_coverage, num_matches, raw_thresholds, coverage_thresholds, act_conf, act_conf_coverage, alt_mut, alt_mut_cover, nontriv_flags = hr.hypothesis_recovery(ref_matrix, sample_vector, ksize, significance=significance, ani_thresh=ani_thresh, min_coverage=min_coverage)
    
    recov_org_data['nontrivial_overlap'] = nontriv_flags
    recov_org_data['in_sample_est'] = is_present
    recov_org_data['num_exclusive_kmers'] = nu
    recov_org_data['num_exclusive_kmers_with_coverage'] = nu_coverage
    recov_org_data['num_matches'] = num_matches
    recov_org_data['acceptance_threshold_wo_coverage'] = raw_thresholds
    recov_org_data['acceptance_threshold_with_coverage'] = coverage_thresholds
    recov_org_data['actual_confidence_wo_coverage'] = act_conf
    recov_org_data['actual_confidence_w_coverage'] = act_conf_coverage
    recov_org_data['p_vals'] = p_vals
    recov_org_data['alt_confidence_mut_rate'] = alt_mut
    recov_org_data['alt_confidence_mut_rate_coverage'] = alt_mut_cover
   
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
    ref_matrix = args.ref_matrix
    sample_file = args.sample_file
    ksize = args.ksize
    ani_thresh = args.ani_thresh
    significance = args.significance
    min_coverage = args.min_coverage
    outfile = args.outfile

    """
      Runs linear program for unknown estimation off of files generated by ref_matrix.py and creates human-readable results file.
      :param matrix_file: location of ref_matrix_processed.npz file (A matrix)
      :param sample_file: location of sample.sig file (y vector)
      :param ksize: kmer size
      :param ani_thresh: ANI cutoff for species equivalence
      :param significance: Minimum probability of individual true negative.
      :param num_kmers_quantile: quantile for determining representative number of kmers in sketch to be used in calculation of p-value.
      :param output_filename: destination for results file; if blank, no file will be written
      :param w: false positive weight. Optional; if set, overrides significance for method 'lp'.
      :return: pandas dataframe containing recovered abundances and metadata.
      """

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
    sample_vector, num_hash_diff_unique, num_hash_diff_total = sample_vector_from_signature(sample_sig, hash_to_idx)
    sample_scale = sample_sig.minhash.scaled
    num_sample_kmers = utils.get_num_kmers(sample_sig, scale=False)
    num_unique_sample_kmers = len(list(sample_sig.minhash.hashes))

    recov_org_data = recover_abundance_data_hyp(
        reference_matrix,
        sample_vector,
        organism_data,
        ksize,
        ani_thresh,
        significance,
        min_coverage,
        num_sample_kmers,
        num_unique_sample_kmers,
        sample_scale
    )

    if outfile:
        recov_org_data.to_csv(outfile)
