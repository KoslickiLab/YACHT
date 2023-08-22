# This will be a flat file of the various functions in various files that are called by the main scripts
# To be used by line_profiler (kernprof -l --view for_profiling.py)
from line_profiler_pycharm import profile
import os
import numpy as np
import pandas as pd
import sys
# add one level up to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
import srcs.hypothesis_recovery_src as hr
from scipy.sparse import load_npz
import argparse
import srcs.utils as utils
import warnings
import os
warnings.filterwarnings("ignore")


script_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))  # currently one level
# above ./tests
data_dir = "../testdata"
out_prefix = "unittest_"
full_out_prefix = os.path.join(data_dir, out_prefix)
abundance_file = full_out_prefix + "recovered_abundance.csv"
reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
sample_sketches = os.path.join(data_dir, "sample.sig")

#############################################################
# py
import numpy as np
import pickle
import sourmash


def load_hashes(filename):
    """
    Helper function that loads the hash_to_col_idx.csv file and returns a dictionary mapping hashes to indices in the
    training dictionary. filename should point to a CSV file with two columns: hash, col_idx.
    :param filename: string (location of the hash_to_col_idx.pkl file)
    :return: dictionary mapping hashes to indicies
    """
    with open(filename, mode='rb') as fid:
        hashes = pickle.load(fid)
    return hashes

@profile
def load_signature_with_ksize(filename, ksize):
    """
    Helper function that loads the signature for a given kmer size from the provided signature file. Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file)
    :param ksize: kmer size
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    return list(sourmash.load_file_as_signatures(filename, ksize=ksize))[0]

@profile
def signatures_mismatch_ksize(signatures, ksize):
    """
    Helper function that checks if any of the signatures in a list have a different kmer size than the given kmer size.
    :param signatures: sourmash signatures
    :param ksize: kmer size
    :return: False (if all signatures have the same kmer size) or the first signature with a different kmer size
    """
    sketch_with_ksize_exists = False
    for sig in signatures:
        if sig.minhash.ksize == ksize:
            sketch_with_ksize_exists = True
            return sketch_with_ksize_exists
    return sketch_with_ksize_exists

@profile
def get_num_kmers(signature, scale=True):
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param signature: sourmash signature
    :return: int (estimated total number of kmers)
    """
    # Abundances may not have been kept, in which case, just use 1
    if signature.minhash.mean_abundance:
        num_kmers = signature.minhash.mean_abundance * len(signature.minhash.hashes)
    else:
        num_kmers = len(signature.minhash.hashes)
    if scale:
        num_kmers *= signature.minhash.scaled
    return np.round(num_kmers)


##################################################
# run_YACHT.py
@profile
def run_YACHT():
    #cmd = f"python {os.path.join(script_dir, 'run_YACHT.py')} --ref_matrix {full_out_prefix}ref_matrix_processed.npz
    # --sample_file " f"{sample_sketches} --outfile {abundance_file} --ksize 31"
    # parse the arguments
    ref_matrix = f"{full_out_prefix}ref_matrix_processed.npz"
    sample_file = sample_sketches
    ksize = 31
    ani_thresh = 0.95
    significance = 0.99
    min_coverage = 1
    outfile = abundance_file

    # check that ksize is an integer
    if not isinstance(ksize, int):
        raise ValueError('ksize must be an integer.')
    # check if min_coverage is between 0 and 1
    if min_coverage < 0 or min_coverage > 1:
        raise ValueError('min_coverage must be between 0 and 1.')

    # Get the training data names
    prefix = ref_matrix.split('ref_matrix_processed.npz')[0]
    hash_to_idx_file = prefix + 'hash_to_col_idx.pkl'
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
    hash_to_idx = load_hashes(hash_to_idx_file)
    organism_data = pd.read_csv(processed_org_file)

    # get the sample y vector (indexed by hash/k-mer, with entry = number of times k-mer appears in sample)
    sample_sig = load_signature_with_ksize(sample_file, ksize)
    # total number of hashes in the training dictionary
    hash_to_idx_keys = list(hash_to_idx.keys())
    K = len(hash_to_idx_keys)
    # initialize the sample vector
    sample_vector = np.zeros(K)
    # get the hashes in the signature (it's for a single sample)
    sample_hashes = sample_sig.minhash.hashes
    # get the hashes that are in both the sample and the training dictionary
    sample_intersect_training_hashes = np.intersect1d(sample_hashes, hash_to_idx_keys, assume_unique=True)
    for sh in sample_intersect_training_hashes:
        idx = hash_to_idx[sh]
        sample_vector[idx] = sample_hashes[sh]

    # get the number of kmers in the sample from the scaled sketch
    sample_scale = sample_sig.minhash.scaled
    #num_sample_kmers = get_num_kmers(sample_sig, scale=False)
    # get the number of unique kmers in the sample
    num_unique_sample_kmers = len(sample_hashes)

    # prep the output data structure, copying over the organism data
    recov_org_data = organism_data.copy()
    #recov_org_data['num_total_kmers_in_sample_sketch'] = num_sample_kmers
    recov_org_data['num_exclusive_kmers_in_sample_sketch'] = num_unique_sample_kmers
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

    # get all the column names of hyp_recovery_df
    hyp_recovery_df_cols = list(hyp_recovery_df.columns)
    # for each of the columns, add it to the recov_org_data
    for col in hyp_recovery_df_cols:
        recov_org_data[col] = hyp_recovery_df[col]

    # TODO: remove the rows that have no overlap with the sample

    # save the results
    recov_org_data.to_csv(outfile)

run_YACHT()