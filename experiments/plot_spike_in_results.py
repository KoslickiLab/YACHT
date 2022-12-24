#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
# import the results
#file = "spike_in_ave_detect.csv"
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="The input CSV file")
args = parser.parse_args()
file = args.input_file
column_header = None
row_header = None
column_indices = []
row_indices = []
values = []
with open(file, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
            column_header = line.strip().split(',')[1]
            row_header = line.strip().split(',')[0]
        elif i == 1:
            column_indices = line.strip().split(',')[2:]
        else:
            row_indices.append(line.strip().split(',')[1])
            values.append(line.strip().split(',')[2:])
# convert the column indices to floats
column_indices = [float(x) for x in column_indices]
# convert the row indices to floats
row_indices = [float(x) for x in row_indices]
# convert the values to floats
values = np.array([[float(x) for x in row] for row in values])

plt.figure()
# make a dataframe
df = pd.DataFrame(values, index=row_indices, columns=column_indices)
# make a heatmap
sns.heatmap(df, cmap="YlGnBu", xticklabels=column_indices, yticklabels=row_indices)
# add a column header
plt.xlabel("--min_coverage parameter")
# add a row header
plt.ylabel("Coverage in sample")
# add the values in the squares
for i in range(len(row_indices)):
    for j in range(len(column_indices)):
        plt.text(j+0.5, i+0.5, values[i, j], horizontalalignment='center', verticalalignment='center')
plt.title("Percentage of times spike-in genome was detected")
# pad the figure with more whitespace
plt.tight_layout()
plt.show()
# save the figure
plt.savefig(f"{file}.png", dpi=300)

