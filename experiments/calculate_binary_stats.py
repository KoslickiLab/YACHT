# This script will compare the TP, FP, and FN of the binary classification of sourmash and our approach
# using the ground truth as baseline
import os
from os.path import exists
import argparse


def main():
    # single argument: directory containing the results of the simulation
    parser = argparse.ArgumentParser(description='Calculate the binary stats for both sourmash and our approach')
    parser.add_argument('-d', '--dir', help='Directory containing the unknown and known signatures', required=True)
    args = parser.parse_args()
    sims_dir = args.dir

    #sims_dir = '.'
    actual_unknowns_file = os.path.join(sims_dir, 'actual_unknowns.txt')
    gather_results_file = os.path.join(sims_dir, 'gather_results.csv')
    our_results_file = os.path.join(sims_dir, 'EU_results_default.csv')
    names_in_sim_file = os.path.join(sims_dir, 'all_names_in_sim.txt')
    reference_names_file = os.path.join(sims_dir, 'known_names.txt')
    if not exists(actual_unknowns_file):
        raise FileNotFoundError(f'Could not find {actual_unknowns_file}. Invoke this script from where the file resides.')
    if not exists(gather_results_file):
        raise FileNotFoundError(f'Could not find {gather_results_file}. Invoke this script from where the file resides.')
    if not exists(our_results_file):
        raise FileNotFoundError(f'Could not find {our_results_file}. Invoke this script from where the file resides.')

    # Get the actual unknowns
    actual_unknowns = []
    with open(actual_unknowns_file, 'r') as f:
        for line in f.readlines():
            actual_unknowns.append(line.strip())
    actual_unknowns = set(actual_unknowns)

    # get the sourmash results
    sourmash_name_loc = 9
    sourmash_results = []
    with open(gather_results_file, 'r') as f:
        # skip the header
        f.readline()
        for line in f.readlines():
            name = line.strip().split(',')[sourmash_name_loc]
            name = name.split('.')[0]
            sourmash_results.append(name)
    sourmash_results = set(sourmash_results)

    # get our results
    our_results = []
    our_name_loc = 1
    rel_abund_loc = 14
    with open(our_results_file, 'r') as f:
        # skip the header
        line = f.readline()
        print(f"Reading our results, using header: {line.strip().split(',')[rel_abund_loc]}")
        for line in f.readlines():
            name = line.strip().split(',')[our_name_loc]
            name = name.split('.')[0]
            rel_abund = float(line.strip().split(',')[rel_abund_loc])
            if rel_abund > 0:
                our_results.append(name)
    our_results = set(our_results)

    # Import the thought to be knowns
    ground_truths = []
    with open(names_in_sim_file, 'r') as f:
        for line in f.readlines():
            line = line.strip().split('.')[0]
            ground_truths.append(line)
    ground_truths = set(ground_truths)

    # remove the unknowns from the ground truths
    ground_truths = ground_truths - actual_unknowns

    # get the names of the organisms in the reference database
    reference_names = []
    with open(reference_names_file, 'r') as f:
        for line in f.readlines():
            reference_names.append(line.strip().split('.')[0])
    reference_names = set(reference_names)

    # get the true negatives
    true_negatives = reference_names - ground_truths

    # calculate TP, FP, and FN for sourmash and our approach
    TP_sourmash = len(ground_truths.intersection(sourmash_results))
    FP_sourmash = len(sourmash_results - ground_truths)
    FN_sourmash = len(ground_truths - sourmash_results)
    TN_sourmash = len(true_negatives) - TP_sourmash - FP_sourmash

    TP_our = len(ground_truths.intersection(our_results))
    FP_our = len(our_results - ground_truths)
    FN_our = len(ground_truths - our_results)
    TN_our = len(true_negatives) - TP_our - FP_our

    sourmash_intersect_our = sourmash_results.intersection(our_results)
    TP_intersect = len(ground_truths.intersection(sourmash_intersect_our))
    FP_intersect = len(sourmash_intersect_our - ground_truths)
    FN_intersect = len(ground_truths - sourmash_intersect_our)
    TN_intersect = len(true_negatives) - TP_intersect - FP_intersect

    # print the results
    print("Method, TP, FP, FN, TN")
    print(f"sourmash, {TP_sourmash}, {FP_sourmash}, {FN_sourmash}, {TN_sourmash}")
    print(f"ours, {TP_our}, {FP_our}, {FN_our}, {TN_our}")
    print(f"intersect, {TP_intersect}, {FP_intersect}, {FN_intersect}, {TN_intersect}")

    # save the results to a file
    with open(os.path.join(sims_dir, 'binary_stats.csv'), 'w') as f:
        f.write("Method,TP,FP,FN,TN\n")
        f.write(f"sourmash,{TP_sourmash},{FP_sourmash},{FN_sourmash},{TN_sourmash}\n")
        f.write(f"ours,{TP_our},{FP_our},{FN_our},{TN_our}\n")
        f.write(f"intersect,{TP_intersect},{FP_intersect},{FN_intersect},{TN_intersect}\n")


if __name__ == '__main__':
    main()
