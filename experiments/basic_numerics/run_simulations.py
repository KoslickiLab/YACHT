import numpy as np
import random
from scipy.sparse import csc_matrix, save_npz
from scipy.stats import binom
from sys import maxsize
import csv
import argparse
import single_sim
import os
from memory_profiler import memory_usage

def run_simulations(
    num_sims,
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
    recovery_method='lp',
    seed=None,
):
    if recovery_method not in {'lp','h'}:
        raise ValueError('Unsupported recovery_method. Currently supported inputs are \'lp\' (linear program) and \'h\' (hypothesis testing)')
    
    if seed is not None:
        np.random.seed(seed)
        sim_seeds = np.random.randint(2**32-1, size=num_sims)
    
    simulation_results = []
    for i in range(num_sims):
        args = (ksize, num_kmers, num_genomes, s_known, s_unknown)
        kwargs = {'mut_thresh':mut_thresh, 'relation_thresh':relation_thresh, 'p_val':p_val, 'coverage':coverage, 'mut_range':mut_range, 'abundance_range':abundance_range,'recovery_method':recovery_method,'seed':sim_seeds[i] if seed is not None else None}
        
        memlist, curr_result = memory_usage(proc=(single_sim.single_sim, args, kwargs), interval=0.01, retval=True)
        maxmem = max(memlist)
        curr_stats = curr_result[0] + [maxmem]
        simulation_results.append((curr_stats, curr_result[1:]))
        
    return simulation_results
    
        
def write_args(filename, args):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header =['num_sims','ksize','num_kmers','num_genomes','s_known','s_unknown','mut_thresh','relation_thresh','p_val','min_coverage','mut_range','abundance_range','recovery_method','seed']
        writer.writerow(header)
        writer.writerow([args.num_sims, args.ksize, args.num_kmers, args.num_genomes, args.s_known, args.s_unknown, args.mut_thresh, args.relation_thresh, args.p_val, args.min_coverage, args.mut_range, args.abundance_range, args.recovery_method, args.seed])
        
        
def write_results(filename, results):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header =['false_positives', 'false_negatives', 'high_fp_mut', 'low_fn_mut', 'n_severe_fps', 'simulation_time', 'recovery_time', 'max_memory_usage']
        writer.writerow(header)
        for res in results:
            writer.writerow(res[0])
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script runs simulations of metagenomic organism detection using a toy model.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num_sims', type=int, help='number of simulations to run')
    parser.add_argument('--output_folder', help='destination for output files', required=True)
    parser.add_argument('--ksize', type=int, help='kmer size', required=True)
    parser.add_argument('--num_kmers', type=int, help='Number of kmers in sketch for each organism', required=True)
    parser.add_argument('--num_genomes', type=int, help='Number of genomes in reference', required=True)
    parser.add_argument('--s_known', type=int, help='Number of known organisms in sample', required=True)
    parser.add_argument('--s_unknown', type=int, help='Number of unknown organisms in sample', required=True)
    parser.add_argument('--mut_thresh', type=float, help='mutation cutoff for species equivalence.', required=False, default = 0.05)
    parser.add_argument('--relation_thresh', type=float, help='relation % of consecutive columns in reference.', required=False, default = 0.95)
    parser.add_argument('--p_val', type=float, help='Maximum probability of individual false negative.', required=False, default = 0.01)
    parser.add_argument('--mut_range', type=float, nargs=2, help='Mutations will be chosen according to uniform distribution within this range. Should satisfy lower < mut_thresh < upper.', required=False, default=[0.01,0.09])
    parser.add_argument('--abundance_range', type=int, nargs=2, help='Abundances will be chosen according to uniform distribution within this range.', required=False, default=[10,101])
    parser.add_argument('--min_coverage', type=float, help='p_val will be valid for organisms with at least this minimum coverage. Should be between 0 and 1.', required=False, default = 1)
    parser.add_argument('--recovery_method', help='Method for recovering organisms; choices are \'lp\' for linear program and \'h\' for hypothesis testing.', required=False, default = 'lp')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility.', default=None, required=False)
    args = parser.parse_args()
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)    
    write_args(args.output_folder + '/args.csv', args)  
    
    sim_results = run_simulations(
        args.num_sims,
        args.ksize,
        args.num_kmers,
        args.num_genomes,
        args.s_known,
        args.s_unknown,
        mut_thresh=args.mut_thresh,
        relation_thresh=args.relation_thresh,
        p_val=args.p_val,
        coverage=args.min_coverage,
        mut_range=args.mut_range,
        abundance_range=args.abundance_range,
        recovery_method=args.recovery_method,
        seed=args.seed,
    )
    
    write_results(args.output_folder + '/results.csv', sim_results)