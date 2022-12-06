import numpy as np
import random
from scipy.sparse import csc_matrix, save_npz
from scipy.stats import binom
from sys import maxsize
import csv
import argparse
import single_sim
import os


def run_simulations(
    n_sims,
    ksize,
    num_unique,
    num_shared,
    num_genomes,
    p_shared,
    s_known,
    s_unknown,
    mut_thresh = 0.05,
    p_val = 0.01,
    coverage = 1,
    mut_range = [0.01,0.09],
    abundance_range = [10,101],
    seed=None,
):
    if seed is not None:
        np.random.seed(seed)
        sim_seeds = np.random.randint(2**32-1, size=n_sims)
    
    simulation_results = []
    for i in range(n_sims):
        curr_result = single_sim.single_sim(
            ksize,
            num_unique,
            num_shared,
            num_genomes,
            p_shared,
            s_known,
            s_unknown,
            mut_thresh = mut_thresh,
            p_val = p_val,
            coverage = coverage,
            mut_range = mut_range,
            abundance_range = abundance_range,
            seed=sim_seeds[i] if seed is not None else None,
        )

        simulation_results.append(curr_result)
        
    return simulation_results
    
        
def write_args(filename, args):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header =['n_sims','ksize','num_unique','num_shared','num_genomes','p_shared','s_known','s_unknown','mut_thresh','p_val','min_coverage','mut_range','abundance_range','seed']
        writer.writerow(header)
        writer.writerow([args.n_sims, args.ksize, args.num_unique, args.num_shared, args.num_genomes, args.p_shared, args.s_known, args.s_unknown, args.mut_thresh, args.p_val, args.min_coverage, args.mut_range, args.abundance_range, args.seed])
        
        
def write_results(filename, results):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        header =['false_positives', 'false_negatives', 'high_fp_mut', 'low_fn_mut', 'n_severe_fps']
        writer.writerow(header)
        for res in results:
            writer.writerow(res[0])
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script runs simulations of metagenomic organism detection using a toy model.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_sims', type=int, help='number of simulations to run')
    parser.add_argument('--output_folder', help='destination for output files', required=True)
    parser.add_argument('--ksize', type=int, help='kmer size', required=True)
    parser.add_argument('--num_unique', type=int, help='Number of kmers unique to each organism', required=True)
    parser.add_argument('--num_shared', type=int, help='Number of kmers shared by each organism', required=True)
    parser.add_argument('--num_genomes', type=int, help='Number of genomes in reference', required=True)
    parser.add_argument('--p_shared', type=float, help='Probability each genome has each shared kmer', required=True)
    parser.add_argument('--s_known', type=int, help='Number of known organisms in sample', required=True)
    parser.add_argument('--s_unknown', type=int, help='Number of unknown organisms in sample', required=True)
    parser.add_argument('--mut_thresh', type=float, help='mutation cutoff for species equivalence.', required=False, default = 0.05)
    parser.add_argument('--p_val', type=float, help='Maximum probability of individual false negative.', required=False, default = 0.01)
    parser.add_argument('--mut_range', type=float, nargs=2, help='Mutations will be chosen according to uniform distribution within this range. Should satisfy lower < mut_thresh < upper.', required=False, default=[0.01,0.09])
    parser.add_argument('--abundance_range', type=int, nargs=2, help='Abundances will be chosen according to uniform distribution within this range.', required=False, default=[10,101])
    parser.add_argument('--min_coverage', type=float, help='p_val will be valid for organisms with at least this minimum coverage. Should be between 0 and 1.', required=False, default = 1)
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility.', default=None, required=False)
    args = parser.parse_args()
    
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)    
    write_args(args.output_folder + '/args.csv', args)  
    
    sim_results = run_simulations(
        args.n_sims,
        args.ksize,
        args.num_unique,
        args.num_shared,
        args.num_genomes,
        args.p_shared,
        args.s_known,
        args.s_unknown,
        mut_thresh=args.mut_thresh,
        p_val=args.p_val,
        coverage=args.min_coverage,
        mut_range=args.mut_range,
        abundance_range=args.abundance_range,
        seed=args.seed,
    )
    
    write_results(args.output_folder + '/results.csv', sim_results)