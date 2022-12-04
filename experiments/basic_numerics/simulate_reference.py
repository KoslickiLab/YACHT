import numpy as np
import random
from scipy.sparse import csc_matrix, save_npz
from scipy.stats import binom
import csv
import argparse


def simulate_ref(num_unique, num_shared, p_shared, num_genomes, seed = None):
    if seed:
        np.random.seed(seed)
    values = [1]*num_unique*num_genomes
    col_idx = []
    row_idx = []
    for i in range(num_genomes):
        #unique kmers
        col_idx += [i]*num_unique
        start = i*num_unique
        end = (i+1)*num_unique
        row_idx += list(range(start,end))
        #shared kmers
        shared_indicators = np.random.binomial(1,p_shared,num_shared)
        shared = list(num_unique*num_genomes + np.nonzero(shared_indicators)[0])
        col_idx += [i]*len(shared)
        row_idx += shared
        values += [1]*len(shared)
    
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
