import numpy as np
import sourmash
import utils


def sample_vector_from_signature(signature, hash_to_idx):
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
    for sh in sig_hash_overlap:
        idx = hash_to_idx[sh]
        sample_vec[idx] = sig_hashes[sh]
    return sample_vec


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
    return sample_vector_from_signature(sample_sig, hash_to_idx)