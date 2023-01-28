import numpy as np
import simulate_sample as ss
import simulate_reference as sr
import sys
sys.path.append('../..')
import recover_abundance as ra
import hypothesis_recovery as hr
import compute_weight as cw
import time


def single_sim(
    ksize,
    num_kmers,
    num_genomes,
    s_known,
    s_unknown,
    ani_thresh = 0.95,
    relation_thresh = 0.95,
    significance = 0.99,
    coverage = 1,
    ani_range = [0.9,1],
    abundance_range = [10,101],
    recovery_method='lp',
    reference_model='det',
    seed=None,
):
    if recovery_method not in {'lp','h'}:
        raise ValueError('Unsupported recovery_method. Currently supported inputs are \'lp\' (linear program) and \'h\' (hypothesis testing)')
    elif reference_model not in {'det','random'}:
        raise ValueError('Unsupported recovery_method. Currently supported inputs are \'det\' (deterministic) and \'random\'')
    
    #generate data
    sim_start = time.time()
    if reference_model == 'det':
        ref_matrix = sr.simulate_ref_deterministic(num_kmers, num_genomes, ksize, relation_thresh=relation_thresh)
    elif reference_model == 'random':
        ref_matrix = sr.simulate_ref_random(num_kmers, num_genomes, ksize, relation_thresh=relation_thresh,seed=seed)
        
    sample_data = ss.simulate_sample(ref_matrix, ksize, s_known, s_unknown, ani_thresh = ani_thresh, ani_range = ani_range, abundance_range = abundance_range, seed=seed)
    y = sample_data[-1]
    sim_end = time.time()
    
    #recovery
    if recovery_method == 'lp':
        w = cw.compute_weight(ksize, num_kmers, p_val = 1-significance, mut_thresh = 1-ani_thresh, coverage = coverage)[0]
        recov_data = ra.recover_abundance_from_vectors(ref_matrix, y, w)
        
    elif recovery_method == 'h':
        recov_data = hr.hypothesis_recovery(ref_matrix, y, ksize, significance=significance, ani_thresh=ani_thresh, min_coverage=coverage)
        
    recov_end = time.time()
    
    sim_time = sim_end - sim_start
    recov_time = recov_end - sim_end
    
    results = sim_results(recov_data, sample_data, ani_thresh = ani_thresh)
    stats = results[0] + [sim_time, recov_time]
    return stats, results[1:]
    
    
def sim_results(recov_data, sample_data, ani_thresh = 0.95):
    support = sample_data[0]
    support_mut = sample_data[1]
    true_known = support[support_mut <= 1-ani_thresh]
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
        
        
