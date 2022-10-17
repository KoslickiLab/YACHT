import numpy as np
import sourmash
import pandas as pd
import gzip
import csv
import argparse
from sklearn.preprocessing import normalize
from scipy.sparse import csc_matrix, save_npz


# input: list of sm signatures
# output: union of all hashes
def all_hashes(signatures):
    sig_hashes = [list(sig.minhash.hashes.keys()) for sig in signatures]
    ah = sig_hashes[0]
    for sh in sig_hashes[1:]:
        ah = np.union1d(ah, sh)
    return dict(zip(ah, list(range(len(ah)))))


def signatures_to_ref_matrix(signatures):
    values = []
    row_idx = []
    col_idx = []
    N = len(signatures)
    hash_to_idx = all_hashes(signatures)
    for col, sig in enumerate(signatures):
        sig_hashes = sig.minhash.hashes
        for h in sig_hashes:
            idx = hash_to_idx[h]
            values.append(sig_hashes[h])
            row_idx.append(idx)
            col_idx.append(col)
    ref_matrix = csc_matrix((values, (row_idx, col_idx)))
    return ref_matrix, hash_to_idx


def process_reference(raw_ref, corr_thresh, max_thresh):
    uncorr_idx = get_uncorr_idx(raw_ref, corr_thresh)
    uncorr_ref = raw_ref[:, uncorr_idx]
    proc_ref = flatten_reference(uncorr_ref, max_thresh)
    return proc_ref, uncorr_idx


def get_uncorr_idx(ref, corr_thresh):
    norm_ref = normalize(ref, norm='l1', axis=0)
    corrs = norm_ref.transpose() * norm_ref
    uncorr_idx = [0]
    N = norm_ref.shape[1]
    for i in range(1, N):
        corr_flag = False
        for j in uncorr_idx:
            if corrs[i, j] > corr_thresh:
                corr_flag = True
                break
        if not (corr_flag):
            uncorr_idx.append(i)
    return np.array(uncorr_idx).astype(int)


def flatten_reference(ref, max_thresh):
    flat_ref = ref.copy()
    flat_ref[flat_ref > max_thresh] = max_thresh
    return flat_ref


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
    writer.writerow(['organism name', 'original index', 'processed index'])
    for i, idx in enumerate(uncorr_org_idx):
        writer.writerow([signatures[idx].name, idx, i])
    f.close()


# input: filename for sourmash signatures
# output: processed and unprocessed matrix files, hash-to-column-indices file, organism manifest
def reference_matrix_from_signatures(signatures, corr_thresh=None, max_thresh=5, mut_thresh=0.05, out_prefix='',N=None):
    if N is not None:
        signatures = signatures[:N]
    if corr_thresh is None:
        k = signatures[0].minhash.ksize
        corr_thresh = 2 * (1 - mut_thresh) ** k
    
    ref_matrix, hashes = signatures_to_ref_matrix(signatures)
    save_npz(out_prefix + 'ref_matrix_unprocessed.npz', ref_matrix)
    
    processed_ref_matrix, uncorr_org_idx = process_reference(ref_matrix, corr_thresh, max_thresh)
    save_npz(out_prefix + 'ref_matrix_processed.npz', processed_ref_matrix)
    
    write_hashes(out_prefix + 'hash_to_col_idx.csv', hashes)
    write_processed_indices(out_prefix + 'processed_org_idx.csv', signatures, uncorr_org_idx)
    return processed_ref_matrix, ref_matrix, hashes, uncorr_org_idx


def reference_matrix_from_file(filename, corr_thresh=None, max_thresh=5, mut_thresh=0.05, out_prefix='', N=None):
    sigs = list(sourmash.load_file_as_signatures(filename))
    return reference_matrix_from_signatures(
        sigs,
        corr_thresh=corr_thresh,
        max_thresh=max_thresh,
        mut_thresh=mut_thresh,
        out_prefix=out_prefix,
        N=N
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Location of signature file')
    parser.add_argument('--corr_thresh', type=float, default=None, help='Threshold for column similarity')
    parser.add_argument('--max_thresh', type=int, default=5, help='Max value of kmer counts')
    parser.add_argument('--mut_thresh', type=float, default=0.05,
                        help='Mutation rate threshold for unrelated organisms')
    parser.add_argument('--out_prefix', help='Location and prefix for output files')
    parser.add_argument('--N', type=int, help='Number of signatures from file to incorporate into matrix')
    args = parser.parse_args()

    reference_matrix_from_file(args.ref_file, corr_thresh=args.corr_thresh, max_thresh=args.max_thresh,mut_thresh=args.mut_thresh, out_prefix=args.out_prefix, N=args.N)
