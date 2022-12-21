import numpy as np
import random
from scipy.sparse import csc_matrix, save_npz
from scipy.stats import binom
import csv
import argparse


# def simulate_ref_random(n_kmers, n_genomes, k, relation_thresh = 0.95, seed = None):
#     if seed is not None:
#         np.random.seed(seed)
#         random.seed(seed)
    
#     total_kmers = np.round(n_kmers*n_genomes*(1-relation_thresh**k))
#     print(total_kmers)
#     values = [1]*n_kmers*n_genomes
#     col_idx = []
#     row_idx = []
#     for i in range(n_genomes):
        
#         kmers_i = list(np.sort(random.sample(range(int(total_kmers)), n_kmers)))
#         col_idx += [i]*n_kmers
#         row_idx += kmers_i
    
#     ref_matrix = csc_matrix((values, (row_idx, col_idx)))
#     return ref_matrix    


def simulate_ref_random(n_kmers, n_genomes, ksize, relation_thresh = 0.95, seed = None):
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)
    
    total_kmers = np.round(n_kmers*n_genomes*(1-relation_thresh**ksize))
    n_overlap = int(np.round((relation_thresh**ksize)*n_kmers))
    n_non_overlap = n_kmers - n_overlap
    total_kmers_non = total_kmers - n_kmers
    all_kmers = list(range(int(total_kmers)))
    
    values = [1]*n_kmers*n_genomes
    col_idx = [0]*n_kmers
    row_idx = list(np.sort(random.sample(range(int(total_kmers)), n_kmers)))
    prev_row_idx = row_idx
    for i in range(1,n_genomes):
    
        overlap_seed = list(np.sort(random.sample(range(int(n_kmers)), n_overlap)))
        overlap = [prev_row_idx[j] for j in overlap_seed]

        non_overlap_seed = list(np.sort(random.sample(range(int(total_kmers_non)), n_non_overlap)))
        non_prev_row_idx = np.setdiff1d(all_kmers, prev_row_idx)
        non_overlap = list(non_prev_row_idx[non_overlap_seed])
        
        kmers_i = overlap + non_overlap
        col_idx += [i]*n_kmers
        row_idx += kmers_i
        
        prev_row_idx = kmers_i
    
    ref_matrix = csc_matrix((values, (row_idx, col_idx)))
    return ref_matrix    
    

def simulate_ref_deterministic(n_kmers, n_genomes, k, relation_thresh = 0.95):
    
    p = relation_thresh**k
    
    col_idx = []
    row_idx = []
    row_idx += list(range(n_kmers))
    col_idx += [0]*n_kmers
    
    breakpoints = [0,n_kmers]
    prev_sizes = [n_kmers]

    for i in range(1, n_genomes):
        new_sizes = [int((p*size)) for size in prev_sizes]
        new_sizes = new_sizes + [int(n_kmers - np.sum(new_sizes))]
        for (j,b) in enumerate(breakpoints):
            row_idx += list(range(b,b+new_sizes[j]))
            col_idx += [i]*new_sizes[j]
        breakpoints.append(b+new_sizes[i])
        prev_sizes = new_sizes
    
    values = [1]*len(col_idx)
    ref_matrix = csc_matrix((values, (row_idx, col_idx)))
    return ref_matrix    

    
def simulate_proc_file(filename, ref_matrix):
    f = open(filename, 'w', newline='', encoding='utf-8')
    writer = csv.writer(f)
    writer.writerow(['organism_name', 'original_index', 'processed_index', 'num_kmers', 'scale_factor','estimated_total_kmers'])
    scale = 1000
    for i in range(np.shape(ref_matrix)[1]):
        kmers = np.sum(ref_matrix[:,i])
        writer.writerow(['sim_organism_'+str(i),i,i,np.sum(ref_matrix[:,i]),scale,scale*kmers])
    f.close()
    

def generate_sim_ref(output_folder, num_unique, num_shared, p_shared, num_genomes, seed=None):
    ref_matrix = simulate_ref(num_unique, num_shared, p_shared, num_genomes, seed)
    save_npz(output_folder + '/ref_matrix_processed.npz', ref_matrix)
    simulate_proc_file(output_folder + '/processed_org_idx.csv', ref_matrix)
     

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script estimates the abundance of microorganisms from a reference database matrix and metagenomic sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_folder', help='destination for output files', required=True)
    parser.add_argument('--num_uniq_kmers', type=int, help='Number of unique kmers per genome', required=True)
    parser.add_argument('--num_shared_kmers', type=int, help='Number of shared kmers', required=True)
    parser.add_argument('--p_shared', type=float, help='Probability each genome has each shared kmer', required=True)
    parser.add_argument('--num_genomes', type=int, help='Size of kmers used in sketch', required=True)
    parser.add_argument('--seed', type=int, help='Size of kmers used in sketch', default=None, required=False)
    args = parser.parse_args()
    
    generate_sim_ref(args.output_folder, args.num_uniq_kmers, args.num_shared_kmers, args.p_shared, args.num_genomes, args.seed)
