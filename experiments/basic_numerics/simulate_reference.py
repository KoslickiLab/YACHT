import numpy as np
import random
from scipy.sparse import csc_matrix, save_npz
from scipy.stats import binom
import csv
import argparse

        
# def simulate_ref(num_hashes, num_genomes, total_hashes):
#     values = [1]*num_hashes*num_genomes
#     col_idx = [0]*num_hashes*num_genomes
#     row_idx = [0]*num_hashes*num_genomes
#     for i in range(num_genomes):
#         col_idx[i*num_hashes:(i+1)*num_hashes] = ([i]*num_hashes)
#         row_idx[i*num_hashes:(i+1)*num_hashes] = random.sample(range(total_hashes),num_hashes)
#     ref_matrix = csc_matrix((values, (row_idx, col_idx)))
#     return ref_matrix

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

# def max_overlap_cdf(x, num_hashes, num_genomes, total_hashes):
#     """
#     Computes the approximate cdf for size of overlaps between two uniform subsets with num_hashes elements out of total_hashes; approximated by binomial distribution.
#     :param x: input for cdf
#     :param ksize: kmer size
#     :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
#     :param num_genomes: number of different genomes in the sample.
#     :return: probability the above random variable is <= x
#     """
#     p = num_hashes/total_hashes
#     single_overlap_cdf = binom.cdf(x, num_hashes, p)
#     return single_overlap_cdf ** num_genomes


# def find_total_hashes(ksize, num_hashes, num_genomes, quantile = None, mut_thresh = 0.05, tol = 1e-10):
#     """
#     Using binary search, finds the total number of hashes such that the CDF of num_hashes*(1-mut_prob)**ksize equals the specified quantile. Quantile defaults to 1-1/num_genomes.
#     :param ksize: kmer size
#     :param num_hashes: expected number of hashes (kmers) in each sketch, or estimate thereof.
#     :param quantile: target quantile
#     :param mut_thresh: mutation probability threshold for species equivalence.
#     :param num_genomes: number of different genomes in the sample.
#     :param tol: search terminates if sufficiently close solution is found.
#     :return: total_hashes number such that CDF of num_hashes*(1-mut_prob)**ksize equals the specified quantile
#     """
#     quantile = 1 - 1/(num_genomes**2) if quantile is None else quantile
#     lower = 1
#     upper = num_hashes*num_genomes
#     mut_prob = (1-mut_thresh)**ksize
#     x = np.ceil(num_hashes*mut_prob)
#     total_hashes_curr = np.ceil(upper/2)
#     q_curr = -1
#     count = 0
#     while (np.abs(quantile - q_curr) > tol):
#         q_curr = max_overlap_cdf(x, num_hashes, num_genomes, total_hashes_curr)
#         if q_curr > quantile:
#             upper = total_hashes_curr
#             total_hashes_curr = lower + np.ceil((upper - lower)/2)
#         else:
#             lower = total_hashes_curr
#             total_hashes_curr = lower + np.ceil((upper - lower)/2)
#         if upper - lower <= 1:
#             total_hashes_curr = upper
#             break
#     return int(total_hashes_curr)