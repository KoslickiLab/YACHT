#!/usr/bin/env python

# this script will take the outputs of end_to_end_full_ref*.sh and make a boxplot of the binary results
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
# make command line args
import argparse
import matplotlib as mpl
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('-fp', '--file_pattern', type=str, required=True, help='Pattern for the simulation folders. Eg. "sims-uniform*"')
args = parser.parse_args()
file_pattern = args.file_pattern
# find all the simulation folders
sim_folders = [x for x in glob(file_pattern) if os.path.isdir(x)]
# for each folder, get the binary_results.csv file, put it in a list, and then concatenate them all together
binary_results = {}
binary_results['sourmash'] = {}
binary_results['ours'] = {}
binary_results['sourmash']['TP'] = []
binary_results['sourmash']['FP'] = []
binary_results['sourmash']['TN'] = []
binary_results['sourmash']['FN'] = []
binary_results['ours']['TP'] = []
binary_results['ours']['FP'] = []
binary_results['ours']['TN'] = []
binary_results['ours']['FN'] = []
for folder in sim_folders:
    # get the binary results file
    # format looks like:
    #
    # Method,TP,FP,FN,TN
    # sourmash,104,499,15,241
    # ours,104,54,15,686
    # intersect,104,34,15,706
    binary_results_file = os.path.join(folder, 'binary_stats.csv')
    with open(binary_results_file, 'r') as fid:
        header = fid.readline().strip().split(',')
        first_line = fid.readline().strip().split(',')
        second_line = fid.readline().strip().split(',')
        for metric in header[1:]:
            binary_results[first_line[0]][metric].append(float(first_line[header.index(metric)]))
            binary_results[second_line[0]][metric].append(float(second_line[header.index(metric)]))

# import this into a pandas dataframe, one row per simulation, columns of metrics and methods
# first convert this to a list of lists
lol = []
for experiment_number in range(len(binary_results['ours']['TP'])):
    for method in binary_results:
        #for metric in binary_results[method]:
        # if you only want the TP and FP
        for metric in ['TP', 'FP']:
            lol.append([experiment_number, method, metric, binary_results[method][metric][experiment_number]])

binary_results_df = pd.DataFrame(lol, columns=['experiment_number', 'method', 'metric', 'value'])
# make a boxplot, with the x axis being the metric and the y axis being the value, one box for each method

sns.boxplot(x="metric", y="value",
            hue="method", palette=["m", "g"],
            data=binary_results_df)
# save the plot
plt.savefig(f"{file_pattern.replace('*','')}binary_results_boxplot.png")
