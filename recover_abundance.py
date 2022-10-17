import numpy as np
import cvxpy as cp
import csv
import sample_vector as sv
from scipy.sparse import load_npz
import argparse
import utils

#inputs: matrix A, vector y, weight w
#output: estimate vector x and metadata
def recover_abundance_from_vectors(A, y, w):
    K, N = np.shape(A)
    x = cp.Variable(N)
    u = cp.Variable(K)
    v = cp.Variable(K)
    tau = 1/(w+1)
    ones_K = np.ones(K)
    objective = cp.Minimize(
        tau*(ones_K @ u) + (1-tau)*(ones_K @ v)
    )
    constraints = [
        x >= 0,
        u >= 0,
        v >= 0,
        u - v + (A @ x) == y,
    ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver = cp.SCIPY, verbose=False)
    return x.value

def recover_abundance_from_files(matrix_file, sample_file, hash_to_idx_file, processed_organism_file, w, output_filename):
    reference_matrix = load_npz(matrix_file)
    sample_vector = sv.sample_vector_from_files(sample_file, hash_to_idx_file)
    abundance = recover_abundance_from_vectors(reference_matrix, sample_vector, w)
    support = np.nonzero(abundance)
    organisms = utils.load_processed_organisms(processed_organism_file)
    write_abundance_results(abundance, organisms, output_filename)
    return abundance

def write_abundance_results(abundance_vector, organisms, output_filename):
    f = open(output_filename, 'w')
    writer = csv.writer(f)
    writer.writerow(['organism name', 'estimated abundance'])
    for i, org in enumerate(organisms):
        writer.writerow([org, abundance_vector[i]])
    f.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script estimates the abundance of microorganisms from a reference database matrix and metagenomic sample.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Reference database matrix in npz format')
    parser.add_argument('--sample_file', help='Metagenomic sample in .sig format')
    parser.add_argument('--hash_file', help='csv file of hash values in database sketch')
    parser.add_argument('--org_file', help='csv list of organisms in database')
    parser.add_argument('--w', type = float, help='False positive weight')
    parser.add_argument('--outfile', help='csv destination for results')
    args = parser.parse_args()
    
    recover_abundance_from_files(args.ref_file, args.sample_file, args.hash_file, args.org_file, args.w, args.outfile)