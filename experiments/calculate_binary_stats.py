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

# Import the thought to be knowns
ground_truths = []
with open(os.path.join(sims_dir, 'known_names.txt'), 'r') as f:
    for line in f.readlines():
        line = line.strip().split('.')[0]
        ground_truths.append(line)
ground_truths = set(ground_truths)

# remove the unknowns from the ground truths
ground_truths = ground_truths - actual_unknowns

# calculate TP, FP, and FN for sourmash and our approach
TP_sourmash = len(ground_truths.intersection(sourmash_results))
FP_sourmash = len(sourmash_results - ground_truths)
FN_sourmash = len(ground_truths - sourmash_results)

TP_our = len(ground_truths.intersection(our_results))
FP_our = len(our_results - ground_truths)
FN_our = len(ground_truths - our_results)

# print the results
print(f"TP sourmash: {TP_sourmash}")
print(f"FP sourmash: {FP_sourmash}")
print(f"FN sourmash: {FN_sourmash}")

print(f"TP our: {TP_our}")
print(f"FP our: {FP_our}")
print(f"FN our: {FN_our}")
