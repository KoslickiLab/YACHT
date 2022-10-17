import numpy as np
import sourmash
import utils

def sample_vector_from_signature(signature, hash_to_idx):
    K = len(list(hash_to_idx.keys()))
    sample_vec = np.zeros(K)
    sig_hashes = signature.minhash.hashes
    sig_hash_overlap = np.intersect1d(sig_hashes, list(hash_to_idx.keys()))
    for sh in sig_hash_overlap:
        idx = hash_to_idx[sh]
        sample_vec[idx] = sig_hashes[sh]
    return sample_vec

def sample_vector_from_files(sig_filename, hash_filename):
    sample_sig = list(sourmash.load_file_as_signatures(sig_filename))[0]
    hash_to_idx = utils.load_hashes(hash_filename)
    return sample_vector_from_signature(sample_sig, hash_to_idx)