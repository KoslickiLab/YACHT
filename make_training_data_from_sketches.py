#!/usr/bin/env python
import numpy as np
import sourmash
import csv
import argparse
from scipy.sparse import csc_matrix, save_npz
import src.utils as utils
import pickle


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
    # first, get the union of all hashes and map them to indices
    sig_hashes = set(x for sig in signatures for x in list(sig.minhash.hashes.keys()))
    hash_to_idx = dict(zip(sig_hashes, range(len(sig_hashes))))
    sig_values = []  # counts of occurrences of hash in genomes
    for col, sig in enumerate(signatures):
        sig_hashes = sig.minhash.hashes
        for hash in sig_hashes:
            idx = hash_to_idx[hash]
            # sig_values is the count of the hash/kmer in the genome
            sig_values.append(sig_hashes[hash])
            row_idx.append(idx)
            col_idx.append(col)
    ref_matrix = csc_matrix((sig_values, (row_idx, col_idx)))
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
    ref_idx = ref_matrix.nonzero()  # indices of nonzero elements in ref_matrix
    mut_prob = (ani_thresh)**ksize  # probability of mutation

    # binary matrix of nonzero elements in ref_matrix
    binary_ref = csc_matrix(([1]*np.shape(ref_idx[0])[0], ref_idx), dtype=bool)
    # number of hashes in each organism
    sizes = np.array(np.sum(binary_ref, axis=0)).reshape(N)

    # sort organisms by size, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    # binary_ref sorted by size
    binary_ref_bysize = binary_ref[:, bysize]

    # all pairwise intersections
    intersections = binary_ref_bysize.T @ binary_ref_bysize
    # set diagonal to 0, since we don't want to compare an organism to itself
    intersections.setdiag([0]*N)
    
    uncorr_idx_bysize = list(range(N))
    for i in range(N):
        intersections_i = intersections[i, uncorr_idx_bysize]
        # if the largest intersection is greater than the mutation probability, remove the organism
        # (since it's too similar, and is small)
        if np.max(intersections_i) > mut_prob*sizes[bysize[i]]:
            uncorr_idx_bysize.remove(i)

    # uncorr_idx is the indices of the organisms in the reference matrix that are uncorrelated
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
    #f = open(filename, 'w', newline='', encoding='utf-8')
    #writer = csv.writer(f)
    #writer.writerow(['hash', 'index'])
    #for h in hashes:
    #    writer.writerow([h, hashes[h]])
    #f.close()


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
    f = open(filename, 'w', newline='', encoding='utf-8')
    writer = csv.writer(f)
    writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
    for i, idx in enumerate(uncorr_org_idx):
        writer.writerow([signatures[idx].name, idx, i, len(signatures[idx].minhash.hashes), utils.get_num_kmers(signatures[idx], scale=False), signatures[idx].minhash.scaled])
    f.close()


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
    signatures = list(sourmash.load_file_as_signatures(ref_file))

    # downsample if desired
    if N is not None:
        signatures = signatures[:N]

    # check that all signatures have the same ksize as the one provided
    sketch_with_ksize_exists = utils.signatures_mismatch_ksize(signatures, ksize)
    if not sketch_with_ksize_exists:
        raise ValueError('There are no signatures with the provided ksize.')

    # convert signatures to reference matrix (rows are hashes/kmers, columns are organisms)
    # TODO: might not need to save the unprocessed matrix
    ref_matrix, hashes = signatures_to_ref_matrix(signatures)
    save_npz(out_prefix + 'ref_matrix_unprocessed.npz', ref_matrix)

    # remove related organisms: any organisms with ANI > 1-mut_thresh are considered the same organism
    processed_ref_matrix, uncorr_org_idx = get_uncorr_ref(ref_matrix, ksize, ani_thresh)
    save_npz(out_prefix + 'ref_matrix_processed.npz', processed_ref_matrix)

    # write out hash-to-row-indices file
    write_hashes(out_prefix + 'hash_to_col_idx.pkl', hashes)

    # write out organism manifest (original index, processed index, num unique kmers, num total kmers, scale factor)
    write_processed_indices(out_prefix + 'processed_org_idx.csv', signatures, uncorr_org_idx)
    #return processed_ref_matrix, ref_matrix, hashes, uncorr_org_idx