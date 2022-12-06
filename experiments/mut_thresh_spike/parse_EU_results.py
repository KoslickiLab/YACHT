#!/usr/bin/env python
# This script will take the EU results, get the ANI's and make the ANI vs. in/out sample mutation threshold plot
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
# Load the metadata for the spikes
sigs_metadata = pd.read_csv('sigs_descriptions.csv', index_col=0)

# get the ANI results
ani_results = pd.read_csv('in_gtdb_similar_to_EU_not_in_sample_clean.csv', index_col=0)
recovery_data = {}  # will be dict: keys are 1-mutation_threshold, values are tuple (ANI, in/out sample)

results_row_basis = []
with open('../formatted_db_processed_org_idx.csv', 'r') as f:
    # skip the first line
    next(f)
    for line in f.readlines():
        results_row_basis.append(line.strip().split(',')[0].split('.')[0])
basis_to_row = {name: i for i, name in enumerate(results_row_basis)}

# find all the EU results
result_files = glob('EU_on_spikes/*')
print("Parsing results")
i = 0
for result_file in result_files:
    i += 1
    if i % 100 == 0:
        print(f"{i} of {len(result_files)}")
    file_name = os.path.basename(result_file)
    # get the spike number
    spike_num = int(file_name.split('_')[5].split('.')[0])
    mut_thresh = float(file_name.split('_')[4])
    # Find out who the similar organism is
    spiked_organism = sigs_metadata.loc[f'{spike_num}.sig', 'name']
    # Find the ANI of the spiked organism to the reference database
    spiked_organism_ani = ani_results.loc[spiked_organism, 'max_ani']
    # find the organism that is similar to the spiked organism
    recover_val = None
    similar_org = ani_results.loc[spiked_organism, 'max_ani_name']
    # Get a single EU result
    EU_result_file = f'EU_on_spikes/EU_results_mut_thresh_{mut_thresh}_{spike_num}.sig_spike.sig.csv'
    row_to_select = basis_to_row[similar_org]
    with open(EU_result_file, 'r') as f:
        # enumerate the lines
        for i, line in enumerate(f):
            if i == row_to_select + 1:
                #print(f"selected line: {line}")
                #print(f"similar org: {similar_org}")
                recover_val = float(line.strip().split(',')[16])
                #print(f"recover val: {recover_val}")
    # check to see if this organism is in the EU results
    in_sample = False
    if recover_val > 0:
        in_sample = True
    ANI_thresh = 1 - mut_thresh
    if ANI_thresh not in recovery_data:
        recovery_data[ANI_thresh] = []
    recovery_data[ANI_thresh].append((spiked_organism_ani, in_sample))

# Make a scatter plot: the x axis is the ANI of the spiked organism to the reference database, the y axis has two
# values: True or false, depending on if the similar organism was in the sample or not
print("Making plot")
for ANI_thresh in recovery_data:
    print(f"ANI threshold: {ANI_thresh}")
    plt.figure()
    x = []
    y = []
    for spiked_organism_ani, in_sample in recovery_data[ANI_thresh]:
        x.append(spiked_organism_ani)
        y.append(in_sample)
    plt.scatter(x, y, label=f'1-mutation threshold: {ANI_thresh}')
    # change the y axis labels to True or False
    plt.yticks([0, 1], ['False', 'True'])
    plt.xlabel('ANI of spiked organism to reference database')
    plt.ylabel('Similar organism detected')
    # change the x-axis range to be [0.7,1]
    plt.xlim([0.7, 1])
    plt.ylim([-0.1, 1.1])
    # add a vertical dashed red line at the ANI threshold
    plt.axvline(x=ANI_thresh, color='r', linestyle='--')
    #plt.show()
    plt.savefig(f'ANI_vs_in_out_sample_mut_thresh_{ANI_thresh}.png')



