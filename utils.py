import numpy as np
import csv
import sourmash


def load_hashes(filename):
    """
    Helper function that loads the hash_to_col_idx.csv file and returns a dictionary mapping hashes to indices in the
    training dictionary. filename should point to a CSV file with two columns: hash, col_idx.
    :param filename: string (location of the hash_to_col_idx.csv file)
    :return: dictionary mapping hashes to indicies
    """
    with open(filename, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        hashes = {int(rows[0]): int(rows[1]) for rows in reader}
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
    for sig in signatures:
        if sig.minhash.ksize != ksize:
            return sig
    return False

def total_kmers_est(signature):
    """
    Helper function that estimates the total number of kmers in a given sample.
    :param signature: sourmash signature
    :return: int (estimated total number of kmers)
    """
    return np.round(signature.minhash.mean_abundance * signature.minhash.scaled * len(signature.minhash.hashes),4)
    