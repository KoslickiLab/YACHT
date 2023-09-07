import os
import numpy as np
import pickle
import sourmash
from tqdm import tqdm


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

    
def load_signature_with_ksize(filename, ksize):
    """
    Helper function that loads the signature for a given kmer size from the provided signature file. Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file)
    :param ksize: kmer size
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    return list(sourmash.load_file_as_signatures(filename, ksize=ksize))[0]


def signatures_mismatch_ksize(signatures, ksize):
    """
    Helper function that checks if any of the signatures in a list have a different kmer size than the given kmer size.
    :param signatures: sourmash signatures
    :param ksize: kmer size
    :return: False (if all signatures have the same kmer size) or True (the first signature with a different kmer size)
    """
    return next(
        (True for sig in signatures if sig.minhash.ksize != ksize),
        False,
    )


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


def check_file_existence(file_path, error_description):
    """
    Helper function that checks if a file exists. If not, raises a ValueError with the given error description.
    :param file_path: string (location of the file)
    :param error_description: string (description of the error)
    :return: None
    """
    if not os.path.exists(file_path):
        raise ValueError(error_description)


def compute_sample_vector(sample_hashes, hash_to_idx):
    """
    Helper function that computes the sample vector for a given sample signature.
    :param sample_hashes: hashes in the sample signature
    :param hash_to_idx: dictionary mapping hashes to indices in the training dictionary
    :return: numpy array (sample vector)
    """
    # total number of hashes in the training dictionary
    hash_to_idx_keys = list(hash_to_idx.keys())
    
    # initialize the sample vector
    sample_vector = np.zeros(len(hash_to_idx_keys))
    
    # get the hashes that are in both the sample and the training dictionary
    sample_intersect_training_hashes = np.intersect1d(sample_hashes, hash_to_idx_keys,  assume_unique=True)
    
    # fill in the sample vector
    for sh in tqdm(sample_intersect_training_hashes):
        sample_vector[hash_to_idx[sh]] = sample_hashes[sh]

    return sample_vector