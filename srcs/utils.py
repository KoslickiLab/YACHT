import os
import numpy as np
import pickle
import sourmash
from tqdm import tqdm, trange
import scipy as sp 
import csv
import zipfile


def load_hashes(filename):
    """
    Helper function that loads the hash_to_col_idx.pkl file and returns a dictionary mapping hashes to indices in the
    training dictionary. filename should point to a CSV file with two columns: hash, col_idx.
    :param filename: string (location of the hash_to_col_idx.pkl file)
    :return: dictionary mapping hashes to indicies
    """
    with open(filename, mode='rb') as fid:
        hashes = pickle.load(fid)
    return hashes

    
def load_signature_with_ksize(filename, ksize):
    """
    Helper function that loads the signature for a given kmer size from the provided signature file.
    Filename should point to a .sig file. Raises exception if given kmer size is not present in the file.
    :param filename: string (location of the signature file)
    :param ksize: kmer size
    :return: sourmash signature
    """
    # Take the first sample signature with the given kmer size
    return list(sourmash.load_file_as_signatures(filename, ksize=ksize))[0]


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

def signatures_to_ref_matrix(signatures, ksize, signature_count):
    """
    Given signature generator, return a sparse matrix with one column per signature and one row per hash
    (union of the hashes)
    :param signatures: sourmash signatures obtained via sourmash.load_file_as_signatures(<file name>)
    :param ksize: kmer size
    :param signature_count: number of signatures in the sourmash signature file
    :return:
    signature_list: list of sourmash signatures
    ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    hash_to_idx: dictionary mapping hash to row index in ref_matrix
    is_mismatch: False (if all signatures have the same kmer size) or True (the first signature with a different kmer size)
    """
    row_idx = []
    col_idx = []
    sig_values = []
    signature_list = []
    is_mismatch = False

    # Use a dictionary to store hash to index mapping
    hash_to_idx = {}
    next_idx = 0  # Next available index for a new hash

    # Iterate over all signatures
    for col, sig in enumerate(tqdm(signatures, total=signature_count)):
        
        # covert to sourmash list
        signature_list.append(sig)
        
        # check that all signatures have the same ksize as the one provided
        if sig.minhash.ksize != ksize:
            is_mismatch = True
            return signature_list, None, None, is_mismatch
        
        sig_hashes = sig.minhash.hashes
        for hash, count in sig_hashes.items():
            # Get the index for this hash， if new hash, and it and increment next_idx
            idx = hash_to_idx.setdefault(hash, next_idx)
            if idx == next_idx:  # New hash was added
                next_idx += 1

            # Append row, col, and value information for creating the sparse matrix
            row_idx.append(idx)
            col_idx.append(col)
            sig_values.append(count)

    # Create the sparse matrix
    ref_matrix = sp.sparse.csc_matrix((sig_values, (row_idx, col_idx)), shape=(next_idx, len(signature_list)))

    return signature_list, ref_matrix, hash_to_idx, is_mismatch

def count_files_in_zip(zip_path):
    """
    Helper function that counts the number of files in a zip file.
    """
    with zipfile.ZipFile(zip_path, 'r') as z:
        return len(z.namelist())

def get_uncorr_ref(ref_matrix, ksize, ani_thresh):
    """
    Given a reference matrix, return a new reference matrix with only uncorrelated organisms
    :param ref_matrix: sparse matrix with one column per signature and one row per hash (union of the hashes)
    :param ksize: int, size of kmer
    :param ani_thresh: threshold for mutation rate, below which we consider two organisms to be correlated/the same
    :return:
    binary_ref: a new reference matrix with only uncorrelated organisms in binary form (discarding counts)
    uncorr_idx: the indices of the organisms in the reference matrix that are uncorrelated/distinct
    """

    N = ref_matrix.shape[1]  # number of organisms
    immut_prob = ani_thresh ** ksize  # probability of immutation

    # Convert the input matrix to a binary matrix
    # since I don't think we are actually using the counts anywhere, we could probably get a performance
    # boost by just using a binary matrix to begin with
    binary_ref = (ref_matrix > 0).astype(int)
    # number of hashes in each organism
    sizes = np.array(np.sum(binary_ref, axis=0)).reshape(N)

    # sort organisms by size  in ascending order, so we keep the largest organism, discard the smallest
    bysize = np.argsort(sizes)
    # binary_ref sorted by size
    binary_ref_bysize = binary_ref[:, bysize]

    # Compute all pairwise intersections
    intersections = binary_ref_bysize.T.dot(binary_ref_bysize)
    # set diagonal to 0, since we don't want to compare an organism to itself
    intersections.setdiag(0)

    uncorr_idx_bysize = np.arange(N)
    immut_prob_sizes = immut_prob * sizes[bysize]

    for i in trange(N):
        # Remove organisms if they are too similar
        # (Note that we remove organism if there is at least one other organism with intersection > immut_prob_sizes[i])
        if np.max(intersections[i, uncorr_idx_bysize]) > immut_prob_sizes[i]:
            uncorr_idx_bysize = np.setdiff1d(uncorr_idx_bysize, i)

    # Sort the remaining indices, uncorr_idx is now the indices of the organisms in the reference matrix that are uncorrelated
    uncorr_idx = np.sort(bysize[uncorr_idx_bysize])

    return binary_ref[:, uncorr_idx], uncorr_idx


def write_hashes(filename, hashes):
    """
    Write a csv file with the following columns: hash, index
    :param filename: output filename
    :param hashes: dictionary mapping hash to index
    :return: None
    """
    with open(filename, 'wb') as fid:
        pickle.dump(hashes, fid)


def write_processed_indices(filename, signatures, uncorr_org_idx):
    """
    Write a csv file with the following columns: organism_name, original_index, processed_index,
    num_unique_kmers_in_genome_sketch, num_total_kmers_in_genome_sketch, genome_scale_factor.
    :param filename: output filename
    :param signatures: sourmash signatures
    :param uncorr_org_idx: the indices of the organisms in the reference matrix that are uncorrelated
    (via get_uncorr_ref)
    :return: None
    """
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_unique_kmers_in_genome_sketch', 'num_total_kmers_in_genome_sketch', 'genome_scale_factor'])
        for i, idx in enumerate(tqdm(uncorr_org_idx)):
            writer.writerow([signatures[idx].name, idx, i, len(signatures[idx].minhash.hashes), get_num_kmers(signatures[idx], scale=False), signatures[idx].minhash.scaled])
