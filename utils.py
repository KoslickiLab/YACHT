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
    return np.round(signature.minhash.mean_abundance * signature.minhash.scaled * len(signature.minhash.hashes),4)


#crude monte carlo method, to be replaced
def estimate_w_temp(k, num_hashes, mut_thresh=0.05, est_n_orgs=1000, p_val=0.99, n_tests=10000):
    prob = (1 - mut_thresh) ** k
    b = []
    for i in range(n_tests):
        b.append(min(np.random.binomial(num_hashes, prob, (est_n_orgs, 1))))
    min_est = np.quantile(b, p_val)
    print(min_est)
    w = min_est / (num_hashes - min_est)
    return w
