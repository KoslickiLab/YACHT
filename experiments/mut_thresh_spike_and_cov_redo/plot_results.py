#!/usr/bin/env python

# NOTE: I don't need this file any more since I'm using the new, md5 approach, so just run check_for_spike_in_recovery_hyp.sh


# This script will take the EU results, get the ANI's and make the ANI vs. in/out sample mutation threshold plot
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
# make command line args
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--results_file', type=str, required=True, help='Results file')
args = parser.parse_args()
results_file = args.results_file

# import the pandas results file
results = pd.read_csv(results_file, index_col=False, sep='@')

# Make a scatter plot: the x axis is the ANI of the spiked organism to the reference database, the y axis has two
# values: True or false, depending on if the similar organism was in the sample or not

spike_coverages = sorted(results['spike_coverage'].unique(), reverse=True)
coverage_thresholds = sorted(results['coverage_threshold'].unique(), reverse=True)
mutation_thresholds = sorted(results['mut_thresh'].unique(), reverse=True)

print(f"spike_coverages: {spike_coverages}")
print(f"coverage_thresholds: {coverage_thresholds}")
print(f"mutation_thresholds: {mutation_thresholds}")

print("Making plot")
for ANI_thresh in mutation_thresholds:
    ANI_thresh = 1 - ANI_thresh
    for spike_coverage in spike_coverages:
        for coverage_threshold in coverage_thresholds:
            print(f"ANI threshold: {ANI_thresh}")
            print(f"Spike coverage: {spike_coverage}")
            print(f"Coverage threshold: {coverage_threshold}")
            x = []
            y = []
            # calculate the percentage of True's/false's above and below the threshold
            ani_below_and_detected = 0
            ani_above_and_not_detected = 0
            # iterate over the rows of the results file
            print(f"spike_coverage: {spike_coverage}, coverage_threshold: {coverage_threshold}, 1 - ANI_thresh: {1 - ANI_thresh}")
            experiment_counter = 0
            for index, row in results.iterrows():
                if np.isclose(row['spike_coverage'], spike_coverage) and np.isclose(row['coverage_threshold'], coverage_threshold) and np.isclose(row['mut_thresh'], 1 - ANI_thresh):
                    experiment_counter += 1
                    if row['rel_ab'] == np.nan:
                        continue
                    spike_ani = row['max_ani']
                    in_sample = False
                    if np.isclose(row['rel_ab'], 1.0):
                        in_sample = True
                    x.append(row['max_ani'])
                    y.append(row['rel_ab'])
                    if spike_ani < ANI_thresh:
                        if in_sample:
                            ani_below_and_detected += 1
                    # false negative only if the spike ANI is above the threshold and the spike coverage was above the coverage threshold
                    if spike_ani >= ANI_thresh:
                        if not in_sample or coverage_threshold > spike_coverage:
                            ani_above_and_not_detected += 1
            plt.figure()
            plt.scatter(x, y, label=f'spiked organism', alpha=0.03)
            # change the y axis labels to True or False
            plt.yticks([0, 1], ['False', 'True'])
            plt.xlabel('ANI of spiked organism to reference database')
            plt.ylabel('Similar organism detected')
            # change the x-axis range to be [0.7,1]
            plt.xlim([0.7, 1])
            plt.ylim([-0.1, 1.1])
            # add a vertical dashed red line at the ANI threshold
            plt.axvline(x=ANI_thresh, color='r', linestyle='--')
            plt.text(ANI_thresh-.009, 0.4, 'ANI threshold', color='r', rotation=90)
            # add a legend describing what the red line is
            plt.legend(loc='upper left')
            # add a title
            plt.title(f'ANI threshold: {ANI_thresh}, Spike coverage: {spike_coverage}, Coverage threshold: {coverage_threshold}')
            #plt.show()
            plt.savefig(f'ANI_vs_in_out_sample_mut_thresh_{ANI_thresh}_cov_thresh_{coverage_threshold}_spike_cov_{spike_coverage}_hyp.png')
            with open(f'ANI_vs_in_out_sample_mut_thresh_{ANI_thresh}_cov_thresh_{coverage_threshold}_spike_cov_{spike_coverage}_FP_and_FN_hyp.txt', 'w') as f:
                FP = ani_below_and_detected / experiment_counter
                FN = ani_above_and_not_detected / experiment_counter
                TP = 1 - FP
                TN = 1 - FN
                f.write("FP,FN,TP,TN\n")
                f.write(f"{FP},{FN},{TP},{TN}\n")




