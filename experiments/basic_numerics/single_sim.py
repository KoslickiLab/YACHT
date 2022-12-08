import numpy as np
import simulate_sample as ss
import simulate_reference as sr
import sys
sys.path.append('../..')
import recover_abundance as ra
import compute_weight as cw
import time


def single_sim(
    ksize,
    num_kmers,
    num_genomes,
    s_known,
    s_unknown,
    mut_thresh = 0.05,
    relation_thresh = 0.95,
    p_val = 0.01,
    coverage = 1,
    mut_range = [0.01,0.09],
    abundance_range = [10,101],
    seed=None,
):
    w = cw.compute_weight(ksize, num_kmers, p_val = p_val, mut_thresh = mut_thresh, coverage = coverage)[0]
    
    sim_start = time.time()
    ref_matrix = sr.simulate_ref_deterministic(num_kmers, num_genomes, ksize, relation_thresh=relation_thresh)
    sample_data = ss.simulate_sample(ref_matrix, ksize, s_known, s_unknown, mut_thresh = mut_thresh, mut_range = mut_range, abundance_range = abundance_range, seed=seed)
    sim_end = time.time()
    
    recov_data = ra.recover_abundance_from_vectors(ref_matrix, sample_data[-1], w)
    recov_end = time.time()
    
    sim_time = sim_end - sim_start
    recov_time = recov_end - sim_end
    
    results = sim_results(recov_data, sample_data, mut_thresh = mut_thresh)
    stats = results[0] + [sim_time, recov_time]
    return stats, results[1:]
    
    
def sim_results(recov_data, sample_data, mut_thresh = 0.05):
    support = sample_data[0]
    support_mut = sample_data[1]
    true_known = support[support_mut <= mut_thresh]
    recov_known = np.nonzero(recov_data[0])[0]
    false_pos = np.setdiff1d(recov_known, true_known)
    false_neg = np.setdiff1d(true_known, recov_known)

    fp_mut = [-1]*np.shape(false_pos)[0]
    #refers to false positives not corresponding to a mutated organism
    severe_fp = []
    for (i, fp) in enumerate(false_pos):
        found = False
        for (j, org) in enumerate(support):
            if fp == org:
                fp_mut[i] = support_mut[j]
                found = True
                break
        if not found:
            severe_fp.append(fp)
        
    fn_mut = [1]*np.shape(false_neg)[0]
    for (i, fn) in enumerate(false_neg):
        for (j, org) in enumerate(support):
            if fn == org:
                fn_mut[i] = support_mut[j]
                break
                
    n_fp = np.shape(false_pos)[0]
    n_fn = np.shape(false_neg)[0]
    high_fp_mut = np.max(fp_mut) if n_fp > 0 else None
    low_fn_mut = np.min(fn_mut) if n_fn > 0 else None
    n_severe_fp = len(severe_fp)
    stats = [n_fp, n_fn, high_fp_mut, low_fn_mut, n_severe_fp]

    return stats, false_pos, false_neg, fp_mut, fn_mut, severe_fp
        
        