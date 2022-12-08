import numpy as np
from scipy.sparse import csc_matrix, save_npz, load_npz
import random


def deletion_matrix(n_rows, nondel_probs):
    n_cols = np.shape(nondel_probs)[0]
    row_idx = []
    col_idx = []
    values = []
    for (i, p) in enumerate(nondel_probs):
        nondel_indicators = np.random.binomial(1, p, n_rows)
        nondel_rows = list(np.nonzero(nondel_indicators)[0])
        row_idx += nondel_rows
        col_idx += [i]*len(nondel_rows)
        values += [1]*len(nondel_rows)
    del_matrix = csc_matrix((values, (row_idx, col_idx)), shape = (n_rows, n_cols))
    return del_matrix


def simulate_sample(ref_matrix, ksize, s_known, s_unknown, mut_thresh = 0.05, mut_range = [0,0.75], abundance_range = [10,101], seed=None):
    if seed:
        np.random.seed(seed)
        random.seed(seed)
    
    rows, cols = np.shape(ref_matrix)
    num_genomes = int(cols)
    s = s_known + s_unknown
    support = np.sort(random.sample(range(num_genomes),s))
    
    known_mut_rts = np.random.uniform(mut_range[0], mut_thresh, s_known) 
    unk_mut_rts = np.random.uniform(mut_thresh, mut_range[1], s_unknown)
    supp_mut_rts = np.hstack([known_mut_rts, unk_mut_rts])
    supp_nomut_probs = (1-supp_mut_rts)**ksize
    
    sample_ref = ref_matrix[:,support]
    n_rows = np.shape(sample_ref)[0]
    del_matrix = deletion_matrix(n_rows, supp_nomut_probs)
    deleted_sample_ref = sample_ref.multiply(del_matrix)
    
    abundance_counts = np.random.randint(abundance_range[0],abundance_range[1],s)
    sample_vector = deleted_sample_ref @ abundance_counts
    
    return support, supp_mut_rts, supp_nomut_probs, abundance_counts, del_matrix, deleted_sample_ref, sample_vector


    