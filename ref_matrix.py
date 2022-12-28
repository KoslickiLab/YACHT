import numpy as np
import sourmash
import pandas as pd
import gzip
import csv
import argparse
import sklearn.preprocessing as sklp
from scipy.sparse import csc_matrix, save_npz
import utils


# input: list of sm signatures
# output: union of all hashes
def all_hashes(signatures):
    sig_hashes = [list(sig.minhash.hashes.keys()) for sig in signatures]
    ah = sig_hashes[0]
    for sh in sig_hashes[1:]:
        ah = np.union1d(ah, sh)
    return dict(zip(ah, list(range(len(ah)))))


def signatures_to_ref_matrix(signatures):
    row_idx = []
    col_idx = []
    N = len(signatures)
    hash_to_idx = all_hashes(signatures)
    sig_values = []
    for col, sig in enumerate(signatures):
        sig_hashes = sig.minhash.hashes
        for h in sig_hashes:
            idx = hash_to_idx[h]
            sig_values.append(sig_hashes[h])
            row_idx.append(idx)
            col_idx.append(col)
    ref_matrix = csc_matrix((sig_values, (row_idx, col_idx)))
    return ref_matrix, hash_to_idx


def get_uncorr_ref(ref, ksize, mut_thresh):
    N = ref.shape[1]
    ref_idx = ref.nonzero()
    mut_prob = (1-mut_thresh)**ksize
    
    binary_ref = csc_matrix(([1]*np.shape(ref_idx[0])[0], ref_idx), dtype=bool)
    sizes = np.array(np.sum(binary_ref, axis = 0)).reshape(N)
    
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
        writer.writerow([signatures[idx].name, idx, i, len(signatures[idx].minhash.hashes), utils.get_num_kmers(signatures[idx], scale = False), signatures[idx].minhash.scaled])
    f.close()


# input: filename for sourmash signatures
# output: processed and unprocessed matrix files, hash-to-column-indices file, organism manifest
def reference_matrix_from_signatures(signatures, ksize, mut_thresh=0.05, out_prefix='',N=None):
    if N is not None:
        signatures = signatures[:N]
        
    mismatch = utils.signatures_mismatch_ksize(signatures, ksize)
    if mismatch:
        raise ValueError(f'Signature for {mismatch.name} has ksize {mismatch.minhash.ksize} that does not match provided ksize {ksize}.')
    
    ref_matrix, hashes = signatures_to_ref_matrix(signatures)
    save_npz(out_prefix + 'ref_matrix_unprocessed.npz', ref_matrix)
    
    processed_ref_matrix, uncorr_org_idx = get_uncorr_ref(ref_matrix, ksize, mut_thresh)
    save_npz(out_prefix + 'ref_matrix_processed.npz', processed_ref_matrix)
    
    write_hashes(out_prefix + 'hash_to_col_idx.csv', hashes)
    write_processed_indices(out_prefix + 'processed_org_idx.csv', signatures, uncorr_org_idx)
    return processed_ref_matrix, ref_matrix, hashes, uncorr_org_idx


def reference_matrix_from_file(filename, ksize, mut_thresh=0.05, out_prefix='', N=None):
    sigs = list(sourmash.load_file_as_signatures(filename))
    return reference_matrix_from_signatures(
        sigs,
        ksize,
        mut_thresh=mut_thresh,
        out_prefix=out_prefix,
        N=N
    )


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

    reference_matrix_from_file(args.ref_file, args.ksize, mut_thresh=args.mut_thresh, out_prefix=args.out_prefix, N=args.N)
