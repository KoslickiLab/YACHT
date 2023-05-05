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
    #with open(filename, mode='r') as infile:
    #    next(infile)
    #    reader = csv.reader(infile)
    #    hashes = {int(rows[0]): int(rows[1]) for rows in reader}
    return hashes

    
def load_signature_with_ksize(filename, ksize):
    """
    Helper function that loads the signature for a given kmer size from the provided signature file. Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file)
    :param ksize: kmer size
    :return: sourmash signature
    """
    sketches = list(sourmash.load_file_as_signatures(filename))
    for sig in sketches:
        if sig.minhash.ksize == ksize:
            return sig
    raise ValueError(f'File {filename} does not contain sketch for ksize = {ksize}.')


def signatures_mismatch_ksize(signatures, ksize):
    """
    Helper function that checks if any of the signatures in a list have a different kmer size than the given kmer size.
    :param signatures: sourmash signatures
    :param ksize: kmer size
    :return: False (if all signatures have the same kmer size) or the first signature with a different kmer size
    """
    for sig in signatures:
        if sig.minhash.ksize != ksize:
            return sig
    return False


def get_num_kmers(signature, scale=True):
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param signature: sourmash signature
    :return: int (estimated total number of kmers)
    """
    num_kmers = signature.minhash.mean_abundance * len(signature.minhash.hashes)
    if scale:
        num_kmers *= signature.minhash.scaled
    return np.round(num_kmers)
    