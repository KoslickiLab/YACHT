#!/usr/bin/env python
import os, sys
import numpy as np
import sourmash
import csv
import argparse
from scipy.sparse import csc_matrix, save_npz
import srcs.utils as utils
import pickle
from tqdm import tqdm
from loguru import logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO");

def signatures_to_ref_matrix(signatures):
    """
    Given a list of signatures, return a sparse matrix with one column per signature and one row per hash
    (union of the hashes)
    :param signatures: list of signatures obtained via sourmash.load_file_as_signatures(<file name>)
    :return:
    ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    hash_to_idx: dictionary mapping hash to row index in ref_matrix
    """
    row_idx = []
    col_idx = []
    sig_values = []
    
    # Use a dictionary to store hash to index mapping
    hash_to_idx = {}
    next_idx = 0  # Next available index for a new hash
    
    # Iterate over all signatures
    for col, sig in enumerate(tqdm(signatures)):
        sig_hashes = sig.minhash.hashes
        for hash, count in sig_hashes.items():
            # Get the index for this hash， if new hash, and it and increment next_idx
            idx = hash_to_idx.setdefault(hash, next_idx)
            if idx == next_idx:  # New hash was added
                next_idx += 1

            # Append row, col, and value information for creating the sparse matrix
            row_idx.append(idx)
            col_idx.append(col)
            sig_values.append(count)

    ref_matrix = csc_matrix((sig_values, (row_idx, col_idx)), shape=(next_idx, len(signatures)))

    return ref_matrix, hash_to_idx


def get_uncorr_ref(ref_matrix, ksize, ani_thresh):
    """
    Given a reference matrix, return a new reference matrix with only uncorrelated organisms
    :param ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    :param ksize: int, size of kmer
    :param ani_thresh: threshold for mutation rate, below which we consider two organisms to be correlated/the same
    :return:
    binary_ref: a new reference matrix with only uncorrelated organisms in binary form (discarding counts)
    uncorr_idx: the indices of the organisms in the reference matrix that are uncorrelated/distinct
    """

    N = ref_matrix.shape[1]  # number of organisms
    immut_prob = ani_thresh ** ksize  # probability of immutation
    
    # Convert the input matrix to a binary matrix
    # since I don't think we are actually using the counts anywhere, we could probably get a performance
    # boost by just using a binary matrix to begin with
    binary_ref = (ref_matrix > 0).astype(int)
    # number of hashes in each organism
    sizes = np.array(np.sum(binary_ref, axis=0)).reshape(N)
    
    # sort organisms by size  in ascending order, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    # binary_ref sorted by size
    binary_ref_bysize = binary_ref[:, bysize]
    
    # Compute all pairwise intersections
    intersections = binary_ref_bysize.T.dot(binary_ref_bysize)
    # set diagonal to 0, since we don't want to compare an organism to itself
    intersections.setdiag(0)
    
    uncorr_idx_bysize = np.arange(N)
    immut_prob_sizes = immut_prob * sizes[bysize]
    
    for i in range(N):        
        # Remove organisms if they are too similar 
        # (Note that we remove organism if there is at least one other organism with intersection > immut_prob_sizes[i])
        if np.max(intersections[i, uncorr_idx_bysize]) > immut_prob_sizes[i]:
            uncorr_idx_bysize = np.setdiff1d(uncorr_idx_bysize, i)

    # Sort the remaining indices, uncorr_idx is now the indices of the organisms in the reference matrix that are uncorrelated
    uncorr_idx = np.sort(bysize[uncorr_idx_bysize])

    return binary_ref[:, uncorr_idx], uncorr_idx


def write_hashes(filename, hashes):
    """
    Write a csv file with the following columns: hash, index
    :param filename: output filename
    :param hashes: dictionary mapping hash to index
    :return: None
    """
    with open(filename, 'wb') as fid:
        pickle.dump(hashes, fid)


def write_processed_indices(filename, signatures, uncorr_org_idx):
    """
    Write a csv file with the following columns: organism_name, original_index, processed_index,
    num_unique_kmers_in_genome_sketch, num_total_kmers_in_genome_sketch, genome_scale_factor.
    :param filename: output filename
    :param signatures: sourmash signatures
    :param uncorr_org_idx: the indices of the organisms in the reference matrix that are uncorrelated
    (via get_uncorr_ref)
    :return: None
    """
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
        for i, idx in enumerate(tqdm(uncorr_org_idx)):
            writer.writerow([signatures[idx].name, idx, i, len(signatures[idx].minhash.hashes), utils.get_num_kmers(signatures[idx], scale=False), signatures[idx].minhash.scaled])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Location of the Sourmash signature file. '
                                           'This is expected to be in Zipfile format (eg. *.sig.zip)', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers in sketch since Zipfiles '
                                                  'can contain multiple k-mer sizes', required=True)
    parser.add_argument('--ani_thresh', type=float, help='mutation cutoff for species equivalence.',
                        required=False, default=0.95)
    parser.add_argument('--out_prefix', help='Location and prefix for output files.', required=True)
    parser.add_argument('--N', type=int, help='Set N to an integer if you only want to take the first N entries'
                                              'in the ref_file; mainly for testing purposes.', required=False)
    args = parser.parse_args()

    # get the arguments
    ref_file = args.ref_file
    ksize = args.ksize
    ani_thresh = args.ani_thresh
    out_prefix = args.out_prefix
    N = args.N

    # load the signatures
    logger.info(f"Loading signatures from {ref_file}")
    signatures = list(sourmash.load_file_as_signatures(ref_file))

    # downsample if desired
    if N is not None:
        signatures = signatures[:N]

    # check that all signatures have the same ksize as the one provided
    # signatures_mismatch_ksize return False (if all signatures have the same kmer size) or True (the first signature with a different kmer size)
    if utils.signatures_mismatch_ksize(signatures, ksize):
        raise ValueError(f"Not all signatures from sourmash signature file {ref_file} have the given ksize {ksize}")

    # convert signatures to reference matrix (rows are hashes/kmers, columns are organisms)
    logger.info("Converting signatures to reference matrix")
    ref_matrix, hashes = signatures_to_ref_matrix(signatures)

    # remove 'same' organisms: any organisms with ANI > ani_thresh are considered the same organism
    logger.info("Removing 'same' organisms with ANI > ani_thresh")
    processed_ref_matrix, uncorr_org_idx = get_uncorr_ref(ref_matrix, ksize, ani_thresh)
    save_npz(f'{out_prefix}_ref_matrix_processed.npz', processed_ref_matrix)

    # write out hash-to-row-indices file
    logger.info("Writing out hash-to-row-indices file")
    write_hashes(f'{out_prefix}_hash_to_col_idx.pkl', hashes)

    # write out organism manifest (original index, processed index, num unique kmers, num total kmers, scale factor)
    logger.info("Writing out organism manifest")
    write_processed_indices(f'{out_prefix}_processed_org_idx.csv', signatures, uncorr_org_idx)
