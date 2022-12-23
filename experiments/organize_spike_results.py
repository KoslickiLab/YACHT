#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import glob
import os
# two arguments: input directory and output file
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="The directory containing the results")
parser.add_argument("output_file", help="The output file")
args = parser.parse_args()
input_dir = args.input_dir
output_file = args.output_file

file_pattern = 'results_iter_*.csv'
results_list = []
thresholds = [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 0.0009765625]
for file in glob.glob(os.path.join(input_dir, file_pattern)):
    # column headers are: sample_coverage   coverage_threshold   recovered_kmer_abundance
    df = pd.read_csv(file)
    # rows are sample coverages, columns are --min_coverage parameter values
    as_matrix = np.zeros((len(thresholds), len(thresholds)))
    for sample_cov_thresh in thresholds:
        for min_coverage in thresholds:
            # get the row and column indices
            row = thresholds.index(sample_cov_thresh)
            col = thresholds.index(min_coverage)
            # get the value
            value = df.loc[(df['sample_coverage'] == sample_cov_thresh) & (df[' coverage_threshold'] == min_coverage), ' recovered_kmer_abundance'].iloc[0]  # oddly enough, pandas doesn't strip the leading space from the column names
            #print(f"row: {row}")
            #print(f"column: {col}")
            #print(f"value: {value}")
            try:
                as_matrix[row, col] = value
            except:
                continue
    results_list.append(as_matrix)


as_3d_array = np.array(results_list)
# average the number of non-zeros over the 1st dimension (iterates)
average = np.mean(as_3d_array>0, axis=0)
# round the averages to 3 decimal places
average = np.round(average, 3)
# add the thresholds to the top and left of the array
average = np.insert(average, 0, thresholds, axis=0)
# make the thresholds the first column
thresholds_with_zero = np.insert(thresholds, 0, 0)
average = np.insert(average, 0, thresholds_with_zero, axis=1)
# convert this to an array of strings
average = average.astype(str)
column_names = ["min_coverage"] * (len(thresholds)+1)
average = np.insert(average, 0, column_names, axis=0)
row_names = ["sample_coverage"] * (len(thresholds)+2)
average = np.insert(average, 0, row_names, axis=1)
np.savetxt(output_file, average, delimiter=",", fmt="%s")

