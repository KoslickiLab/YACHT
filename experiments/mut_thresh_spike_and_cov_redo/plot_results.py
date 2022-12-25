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
results = pd.read_csv(results_file, index_col=0)

# Make a scatter plot: the x axis is the ANI of the spiked organism to the reference database, the y axis has two
# values: True or false, depending on if the similar organism was in the sample or not
print("Making plot")
for ANI_thresh in recovery_data:
    print(f"ANI threshold: {ANI_thresh}")
    plt.figure()
    x = []
    y = []
    # calculate the percentage of True's/falses above and below the threshold
    ani_below_and_detected = 0
    ani_above_and_not_detected = 0
    for spiked_organism_ani, in_sample in recovery_data[ANI_thresh]:
        x.append(spiked_organism_ani)
        y.append(in_sample)
        if spiked_organism_ani < ANI_thresh:
            if in_sample:
                ani_below_and_detected += 1
        if spiked_organism_ani >= ANI_thresh:
            if not in_sample:
                ani_above_and_not_detected += 1
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
    #plt.show()
    plt.savefig(f'ANI_vs_in_out_sample_mut_thresh_{ANI_thresh}_hyp.png')
    with open(f'ANI_vs_in_out_sample_mut_thresh_{ANI_thresh}_FP_and_FN_hyp.txt', 'w') as f:
        FP = ani_below_and_detected / len(recovery_data[ANI_thresh])
        FN = ani_above_and_not_detected / len(recovery_data[ANI_thresh])
        TP = 1 - FP
        TN = 1 - FN
        f.write("FP,FN,TP,TN\n")
        f.write(f"{FP},{FN},{TP},{TN}\n")




