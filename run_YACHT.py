#!/usr/bin/env python
import os, sys
import numpy as np
import pandas as pd
import srcs.hypothesis_recovery_src as hr
from scipy.sparse import load_npz
import argparse
import srcs.utils as utils
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from loguru import logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO");

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
    if not (0 <= min_coverage <= 1):
        raise ValueError('min_coverage must be between 0 and 1.')

    # Get the training data names
    prefix = ref_matrix.split('ref_matrix_processed.npz')[0]
    hash_to_idx_file = prefix + 'hash_to_col_idx.pkl'
    processed_org_file = prefix + 'processed_org_idx.csv'

    # make sure all these files exist
    utils.check_file_existence(ref_matrix, f'Reference matrix file {ref_matrix} does not exist. Please run ref_matrix.py first.')
    utils.check_file_existence(hash_to_idx_file, f'Hash to index file {hash_to_idx_file} does not exist. Please run ref_matrix.py first.')
    utils.check_file_existence(processed_org_file, f'Processed organism file {processed_org_file} does not exist. Please run ref_matrix.py first.')

    # load the training data
    logger.info('Loading reference matrix, hash to index dictionary, and organism data.')
    reference_matrix = load_npz(ref_matrix)
    hash_to_idx = utils.load_hashes(hash_to_idx_file)
    organism_data = pd.read_csv(processed_org_file)

    logger.info('Loading sample signature.')
    # get the sample y vector (indexed by hash/k-mer, with entry = number of times k-mer appears in sample)
    sample_sig = utils.load_signature_with_ksize(sample_file, ksize)
    
    logger.info('Computing sample vector.')
    # get the hashes in the sample signature (it's for a single sample)
    sample_hashes = sample_sig.minhash.hashes
    sample_vector = utils.compute_sample_vector(sample_hashes, hash_to_idx)

    # get the number of kmers in the sample from the scaled sketch
    num_sample_kmers = utils.get_num_kmers(sample_sig, scale=False)  # TODO: might not save this for time reasons
    # get the number of unique kmers in the sample
    num_unique_sample_kmers = len(sample_hashes)

    # prep the output data structure, copying over the organism data
    recov_org_data = organism_data.copy()
    recov_org_data['num_total_kmers_in_sample_sketch'] = num_sample_kmers  # TODO: might not save this for time reasons
    recov_org_data['num_exclusive_kmers_in_sample_sketch'] = num_unique_sample_kmers
    recov_org_data['sample_scale_factor'] = sample_sig.minhash.scaled
    recov_org_data['min_coverage'] = min_coverage

    # check that the sample scale factor is the same as the genome scale factor for all organisms
    sample_diff_idx = np.where(recov_org_data['sample_scale_factor'].ne(recov_org_data['genome_scale_factor']).to_list())[0].tolist()
    sample_diffs = recov_org_data['organism_name'].iloc[sample_diff_idx]
    if not sample_diffs.empty:
        raise ValueError(f'Sample scale factor does not equal genome scale factor for organism {sample_diffs.iloc[0]} and {len(sample_diffs) - 1} others.')

    # compute hypothesis recovery
    logger.info('Computing hypothesis recovery.')
    hyp_recovery_df, nontriv_flags = hr.hypothesis_recovery(
        reference_matrix, sample_vector, ksize, significance=significance, ani_thresh=ani_thresh, min_coverage=min_coverage)
    
    # Boolean indicating whether genome shares at least one k-mer with sample
    recov_org_data['nontrivial_overlap'] = nontriv_flags
    
    # for each of the columns of hyp_recovery_df, add it to the recov_org_data
    for col in hyp_recovery_df.columns:
        recov_org_data[col] = hyp_recovery_df[col]

    # remove from recov_org_data all those with non-trivial overlap 0
    recov_org_data = recov_org_data[recov_org_data['nontrivial_overlap'] == 1]

    # save the results
    recov_org_data.to_csv(outfile, index=None)
