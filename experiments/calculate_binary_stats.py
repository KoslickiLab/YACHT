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

# get the sourmash results
sourmash_name_loc = 4
sourmash_results = []
with open(os.path.join(sims_dir, 'gather_results.csv'), 'r') as f:
    for line in f.readlines():
        name = line.strip().split(',')[sourmash_name_loc]
        sourmash_results.append(name)

# get our results
our_results = []
our_name_loc = 1
rel_abund_loc = 24
with open(os.path.join(sims_dir, 'EU_results_default.csv'), 'r') as f:
    for line in f.readlines():
        name = line.strip().split(',')[our_name_loc]
        rel_abund = float(line.strip().split(',')[rel_abund_loc])
        if rel_abund >= 0:
            our_results.append(name)


