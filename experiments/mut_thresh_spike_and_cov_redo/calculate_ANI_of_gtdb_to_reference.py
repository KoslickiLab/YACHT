#!/usr/bin/env python
# script to take all the absent organisms and calculate the ANI of them to the reference database
import os
#os.chdir('./experiments/mut_thresh_spike')
import sourmash
import multiprocessing
from itertools import repeat
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Calculate the ANI of the absent organisms to the reference database')
parser.add_argument('--mutation_rate', help='the mutation rate that was used to form the reference database')
parser.add_argument('--gtdb', help='the gtdb .sig file for those not in the sample but similar to the merged reference')
parser.add_argument('--reference_database_full', help='the reference database')
# parse the arguments
args = parser.parse_args()
mutation_rate = args.mutation_rate
reference_name = args.gtdb
EU_training_name = args.reference_database_full

#reference_name = "gtdb-rs207.genomic-reps.dna.k31_not_in_sample_similar_to_merged_reference.sig"
#EU_training_name = "../formatted_db.sig"
# import the reference database
reference_db = sourmash.load_file_as_signatures(reference_name)
ref_sigs = {sketch.name.split('.')[0]: sketch for sketch in reference_db}
# import the absent organisms
EU_training_db = sourmash.load_file_as_signatures(EU_training_name)
EU_training_sigs = {sketch.name.split('.')[0]: sketch for sketch in EU_training_db}

# import the names of the orgs that passed the filter
df = pd.read_csv(f"formatted_db_mut_{mutation_rate}_processed_org_idx.csv", sep=',')
# get the names of the orgs that passed the filter
passed_orgs = df['organism_name'].values
# just take the first part of the accession
passed_orgs = [passed_org.split('.')[0] for passed_org in passed_orgs]


def map_func(ref_sig, EU_training_sketches):
    # map function to calculate the ANI of the absent organisms to the reference database
    max_ani = 0
    max_name = None
    for absent_name, absent_sketch in EU_training_sketches.items():
        ani = absent_sketch.max_containment_ani(ref_sig)
        ani = ani.ani
        if ani and ani > max_ani:
            max_ani = ani
            max_name = absent_name
    return ref_sig.name, ref_sig.md5sum(), max_ani, max_name


def map_star(args):
    # wrapper function to unpack the arguments
    return map_func(*args)

# now we can calculate the ANI of the absent organisms to the reference database in parallel
N = 100
pool = multiprocessing.Pool(processes=N)
results = pool.map(map_star, zip(ref_sigs.values(), repeat(EU_training_sigs)), chunksize=int(len(ref_sigs)/N))
pool.close()
pool.join()
# now we can write the results to a file
with open(f'in_gtdb_similar_to_EU_not_in_sample_mut_{mutation_rate}.tsv', 'w') as f:
    f.write('gtdb_name\tgtdb_md5\tmax_ani\tmax_ani_name\n')
    for ref_name, md5, ani, name in results:
        if ani > 0.7:
            if name in passed_orgs:
                f.write(f'{ref_name}\t{md5}\t{ani}\t{name}\n')


