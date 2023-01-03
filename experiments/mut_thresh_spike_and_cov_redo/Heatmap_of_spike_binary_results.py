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
import matplotlib as mpl
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--results_file', type=str, required=True, help='Results file')
args = parser.parse_args()
results_file = args.results_file

# import the pandas results file
results = pd.read_csv(results_file, index_col=False, sep='@')

# Make a scatter plot: the x axis is the ANI of the spiked organism to the reference database, the y axis has two
# values: True or false, depending on if the similar organism was in the sample or not

spike_coverages = sorted(results['spike_coverage'].unique(), reverse=True)
coverage_thresholds = sorted(results['coverage_threshold'].unique(), reverse=True)
mutation_thresholds = sorted(results['mut_thresh'].unique(), reverse=True)
ANI_thresholds = [1 - x for x in mutation_thresholds]

#print(f"spike_coverages: {spike_coverages}")
#print(f"coverage_thresholds: {coverage_thresholds}")
#print(f"mutation_thresholds: {mutation_thresholds}")

# all the files have the form ANI_vs_in_out_sample_mut_thresh_0.95_cov_thresh_0.01_spike_cov_0.01_FP_and_FN_hyp.txt
# so for each mut_thresh, let's make a 2x2 heatmap of the cov thresh vs the spike cov
# for each mut_thresh, let's make a 2x2 heatmap of the cov thresh vs the spike cov
for ANI_threshold in ANI_thresholds:
    # make an array of the cov_thresh vs spike_cov
    fp_array = np.ones((len(coverage_thresholds), len(spike_coverages)))
    fn_array = np.ones((len(coverage_thresholds), len(spike_coverages)))
    for coverage_threshold_index, coverage_threshold in enumerate(coverage_thresholds):
        for spike_coverage_index, spike_coverage in enumerate(spike_coverages):
            # import the results file
            results_file = f"ANI_vs_in_out_sample_mut_thresh_{ANI_threshold}_cov_thresh_{coverage_threshold}_spike_cov_{spike_coverage}_FP_and_FN_hyp.txt"
            df = pd.read_csv(results_file, sep=',', index_col=False)
            fp_array[coverage_threshold_index, spike_coverage_index] = df['FP'].values[0]
            fn_array[coverage_threshold_index, spike_coverage_index] = df['FN'].values[0]
            print(f"ANI_threshold: {ANI_threshold}\n coverage_threshold: {coverage_threshold}\n spike_coverage: {spike_coverage}\n FP: {df['FP'].values[0]}\n FN: {df['FN'].values[0]}\n")
    # round to 2 decimal places
    fp_array = np.round(fp_array, 3).T
    fn_array = np.round(fn_array, 3).T
    # make the heatmap
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    cmap_type = 'viridis'
    ax[0].imshow(fp_array, cmap=cmap_type, interpolation='nearest', vmin=0, vmax=1, aspect='auto')
    ax[0].set_title(f"False Positives, mut_thresh={ANI_threshold}")
    ax[0].set_yticks(range(len(spike_coverages)))
    ax[0].set_xticks(range(len(coverage_thresholds)))
    # table is flipped horiz and vert, so need to flip the labels
    ax[0].set_yticklabels(spike_coverages)
    ax[0].set_xticklabels(coverage_thresholds)
    ax[0].set_ylabel('spike_coverage')
    ax[0].set_xlabel('coverage_threshold')
    for i in range(len(coverage_thresholds)):
        for j in range(len(spike_coverages)):
            text = ax[0].text(j, i, fp_array[i, j], ha="center", va="center", color="w")
    plt.set_cmap('viridis')
    # add a colorbar
    #im = ax[0].imshow(fp_array, vmin=0, vmax=1)
    #cbar = ax[0].figure.colorbar(im, ax=ax[0])
    #cbar.ax.set_ylabel("False Positives", rotation=-90, va="bottom")
    im = ax[1].imshow(fn_array, cmap=cmap_type, interpolation='nearest', vmin=0, vmax=1, aspect='auto')
    ax[1].set_title(f"False Negatives, mut_thresh={ANI_threshold}")
    ax[1].set_yticks(range(len(spike_coverages)))
    ax[1].set_xticks(range(len(coverage_thresholds)))
    ax[1].set_yticklabels(spike_coverages)
    ax[1].set_xticklabels(coverage_thresholds)
    ax[1].set_ylabel('spike_coverage')
    ax[1].set_xlabel('coverage_threshold')
    for i in range(len(coverage_thresholds)):
        for j in range(len(spike_coverages)):
            text = ax[1].text(j, i, fn_array[i, j], ha="center", va="center", color="w")
    # add a colorbar, scaled to be from 0 to 1
    #im = ax[1].imshow(fn_array, vmin=0, vmax=1)
    cbar = ax[1].figure.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)
    #cbar.ax.set_ylabel("False Negatives", rotation=-90, va="bottom")
    fig.tight_layout()
    # change the color scheme to be more intuitive
    plt.set_cmap(cmap_type)
    plt.savefig(f"heatmap_mut_thresh_{ANI_threshold}.png")
    plt.close()

