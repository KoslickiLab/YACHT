# This script will compare the TP, FP, and FN of the binary classification of sourmash and our approach
# using the ground truth as baseline
import os
import argparse

# Get the actual unknowns
sims_dir = "."
actual_unknowns = []
with open(os.path.join(sims_dir, 'actual_unknowns.txt'), 'r') as f:
    for line in f.readlines():
        actual_unknowns.append(line.strip())
actual_unknowns = set(actual_unknowns)

# get the sourmash results
sourmash_name_loc = 9
sourmash_results = []
with open(os.path.join(sims_dir, 'gather_results.csv'), 'r') as f:
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
rel_abund_loc = 13
with open(os.path.join(sims_dir, 'EU_results_default.csv'), 'r') as f:
    # skip the header
    f.readline()
    for line in f.readlines():
        name = line.strip().split(',')[our_name_loc]
        name = name.split('.')[0]
        rel_abund = float(line.strip().split(',')[rel_abund_loc])
        if rel_abund >= 0:
            our_results.append(name)
our_results = set(our_results)
# calculate TP, FP, and FN for sourmash and our approach
sourmash_TP = len(sourmash_results.intersection(actual_unknowns))
sourmash_FP = len(sourmash_results.difference(actual_unknowns))
sourmash_FN = len(actual_unknowns.difference(sourmash_results))

our_TP = len(our_results.intersection(actual_unknowns))
our_FP = len(our_results.difference(actual_unknowns))
our_FN = len(actual_unknowns.difference(our_results))

# print the results
print(f"Sourmash TP: {sourmash_TP}")
print(f"Sourmash FP: {sourmash_FP}")
print(f"Sourmash FN: {sourmash_FN}")

print(f"Our TP: {our_TP}")
print(f"Our FP: {our_FP}")
print(f"Our FN: {our_FN}")


