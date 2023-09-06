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
    
