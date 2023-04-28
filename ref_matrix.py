import numpy as np
import sourmash
import pandas as pd
import gzip
import csv
import argparse
import sklearn.preprocessing as sklp
from scipy.sparse import csc_matrix, save_npz
import utils


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


def get_uncorr_ref(ref, ksize, mut_thresh):
    N = ref.shape[1]
    ref_idx = ref.nonzero()
    mut_prob = (1-mut_thresh)**ksize
    
    binary_ref = csc_matrix(([1]*np.shape(ref_idx[0])[0], ref_idx), dtype=bool)
    sizes = np.array(np.sum(binary_ref, axis=0)).reshape(N)
    
    bysize = np.argsort(sizes)
    binary_ref_bysize = binary_ref[:,bysize]
    
    intersections = binary_ref_bysize.T @ binary_ref_bysize
    intersections.setdiag([0]*N)
    
    uncorr_idx_bysize = list(range(N))
    for i in range(N):
        intersections_i = intersections[i, uncorr_idx_bysize]
        if np.max(intersections_i) > mut_prob*sizes[bysize[i]]:
            uncorr_idx_bysize.remove(i)
        
    uncorr_idx = np.sort(bysize[uncorr_idx_bysize])
    
    return binary_ref[:,uncorr_idx], uncorr_idx


def write_hashes(filename, hashes):
    f = open(filename, 'w', newline='', encoding='utf-8')
    writer = csv.writer(f)
    writer.writerow(['hash', 'index'])
    for h in hashes:
        writer.writerow([h, hashes[h]])
    f.close()


def write_processed_indices(filename, signatures, uncorr_org_idx):
    f = open(filename, 'w', newline='', encoding='utf-8')
    writer = csv.writer(f)
    writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
    for i, idx in enumerate(uncorr_org_idx):
        writer.writerow([signatures[idx].name, idx, i, len(signatures[idx].minhash.hashes), utils.get_num_kmers(signatures[idx], scale=False), signatures[idx].minhash.scaled])
    f.close()


# input: filename for sourmash signatures
# output: processed and unprocessed matrix files, hash-to-column-indices file, organism manifest


def reference_matrix_from_file(filename, ksize, mut_thresh=0.05, out_prefix='', N=None):
    """
    This function loads a collection of signatures from a file and converts them into a reference matrix.
    :param filename: filename of the signature file (.sig.zip) from sourmash
    :param ksize: what k-size to use (as a zipfile can contain multiple k-sizes)
    :param mut_thresh: 1-ANI threshold, i.e. the threshold for what is considered a different organism (eg. 0.05 means
    organisms with >95% ANI are considered the same)
    :param out_prefix: prefix for output files
    :param N: Integer, if you want to downsample the number of organisms in the reference matrix for testing purposes
    :return:
    processed_ref_matrix: the reference matrix with rows corresponding to hash values and columns corresponding to
    organisms, with only one of any two organisms kept with ANI > 1-mut_thresh
    ref_matrix: the reference matrix with rows corresponding to hash values and columns corresponding to organisms,
    with all organisms kept
    hashes: a dictionary mapping hash values to row indices in the reference matrix
    (i.e. an ordered list of hash values)
    uncorr_org_idx: a list of indices of the organisms in the processed reference matrix that correspond to the
    original organisms in the unprocessed reference matrix. This also includes additional information such as:
    num_unique_kmers_in_genome_sketch,num_total_kmers_in_genome_sketch,genome_scale_factor
    """
    # load the signatures
    signatures = list(sourmash.load_file_as_signatures(filename))
    # downsample if desired
    if N is not None:
        signatures = signatures[:N]

    # check that all signatures have the same ksize as the one provided
    mismatch = utils.signatures_mismatch_ksize(signatures, ksize)
    if mismatch:
        raise ValueError(
            f'Signature for {mismatch.name} has ksize {mismatch.minhash.ksize} that does not match provided ksize {ksize}.')

    # convert signatures to reference matrix (rows are hashes/kmers, columns are organisms)
    # TODO: might not need to save the unprocessed matrix
    ref_matrix, hashes = signatures_to_ref_matrix(signatures)
    save_npz(out_prefix + 'ref_matrix_unprocessed.npz', ref_matrix)

    # remove related organisms: any organisms with ANI > 1-mut_thresh are considered the same organism
    processed_ref_matrix, uncorr_org_idx = get_uncorr_ref(ref_matrix, ksize, mut_thresh)
    save_npz(out_prefix + 'ref_matrix_processed.npz', processed_ref_matrix)

    # write out hash-to-row-indices file
    write_hashes(out_prefix + 'hash_to_col_idx.csv', hashes)
    # write out organism manifest (original index, processed index, num unique kmers, num total kmers, scale factor)
    write_processed_indices(out_prefix + 'processed_org_idx.csv', signatures, uncorr_org_idx)
    return processed_ref_matrix, ref_matrix, hashes, uncorr_org_idx


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Location of signature file', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers in sketch', required=True)
    parser.add_argument('--mut_thresh', type=float, default=0.05,
                        help='Mutation rate threshold for unrelated organisms', required=False)
    parser.add_argument('--out_prefix', help='Location and prefix for output files', required=True)
    parser.add_argument('--N', type=int, help='Number of signatures from file to incorporate into matrix', required=False)
    args = parser.parse_args()

    ref_file = args.ref_file
    ksize = args.ksize
    mut_thresh = args.mut_thresh
    out_prefix = args.out_prefix
    N = args.N
    reference_matrix_from_file(ref_file, ksize, mut_thresh=mut_thresh, out_prefix=out_prefix, N=N)
