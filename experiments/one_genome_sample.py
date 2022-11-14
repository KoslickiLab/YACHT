#!/usr/bin/env python

# This script will get everything in order to run a "one genome" experiment
# prerequisites: Run  sourmash compare --containment --estimate-ani -o formatted_db.sig.pairwise formatted_db.sig
import os
import argparse
import sourmash
import numpy as np
import pandas as pd
import glob
import shutil

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Creates a one genome experiment with the the file names'
                                                 'expected by the single_genome.sh script. Results saved in the output '
                                                 'directory results.txt, and NOT in the main directory too.'
                                                 'By default, there are 26 buckets, with at least 24 pairs in each '
                                                 'bucket. Buckets are 0.2 ANI wide. Smallest bucket goes down to 0.74 '
                                                 'ANI')
    parser.add_argument('-d', '--out_dir', help='Directory containing the unknown and known signatures',
                        required=False, default='one_genome_sample')
    parser.add_argument('-bn', '--bucket_number', help='Which bucket you want to select from. 1-26 are your options.',
                        required=True, type=int)
    parser.add_argument('-pn', '--pair_number', help='Which pair you want to select from. 1 to 24 are safe options.',
                        required=True, type=int)
    args = parser.parse_args()
    output_dir = args.out_dir
    bucket_number = args.bucket_number
    pair_number = args.pair_number

    if not os.path.exists('MANIFEST.csv'):
        raise FileNotFoundError(f"Could not find MANIFEST.csv {os.getcwd()}")

    pairwise_file = 'formatted_db.sig.pairwise'
    if not os.path.exists(pairwise_file):
        raise FileNotFoundError(
            f"Could not find {pairwise_file} in the same directory as where this was run. Please run sourmash compare "
            f"--containment --estimate-ani -o formatted_db.sig.pairwise formatted_db.sig")
    basis_file = f"{pairwise_file}.labels.txt"
    if not os.path.exists(basis_file):
        raise FileNotFoundError(
            f"Could not find {basis_file} in the same directory as where this was run. Please run sourmash compare "
            f"--containment --estimate-ani -o formatted_db.sig.pairwise formatted_db.sig")

    genome_dir = "ref_genomes_3/reference_genomes/"

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)



    # Get all the names
    with open(basis_file, 'r') as f:
        all_names = [line.strip().split('.')[0] for line in f.readlines() if line.strip()]

    index_2_name = {i: name for i, name in enumerate(all_names)}

    # get the pw distances
    pw_distances = np.load(pairwise_file)

    # get the upper triangle of the pw distances
    upper_triangle = np.triu(pw_distances, k=1)

    index_2_value = {}
    for i in range(upper_triangle.shape[0]):
        for j in range(upper_triangle.shape[1]):
            if j > i:
                index_2_value[(i, j)] = upper_triangle[i, j]

    # put these in a data frame, columns are name of the ith corrdinate, name of the jth coordinate, and the distance
    # first, make a dictionary with all this info
    data = {'i': [], 'j': [], 'distance': []}
    for (i, j), distance in index_2_value.items():
        data['i'].append(index_2_name[i])
        data['j'].append(index_2_name[j])
        data['distance'].append(distance)
    # stick it in a data frame
    df = pd.DataFrame(data)
    # sort descending based on distance
    df = df.sort_values(by='distance', ascending=False)

    # make a list of data frames, one for each bucket 1:0.1:91
    df_list = []
    for i in range(26):
        df_list.append(df[(df['distance'] <= 1 - i / 100) & (df['distance'] > 1 - (i + 1) / 100)])

    # check if there are any shared i's in the data frames
    #tops = set(df_list[0]['i'])
    #for i in range(1, 4):
    #    tops = tops.intersection(set(df_list[i]['i']))

    # Looks like GCF_016625365 is a good anchor, it's in the first 4 buckets
    #anchor = 'GCF_016625365'
    reference = df_list[bucket_number].iloc[pair_number]['i']
    simulation = df_list[bucket_number].iloc[pair_number]['j']
    distance = df_list[bucket_number].iloc[pair_number]['distance']

    # write the fake simulated data
    # find the files with the anchor as a prefix
    reference_files = glob.glob(os.path.join(genome_dir, "**", f"{reference}*"), recursive=True)
    to_copy = ''
    for file in reference_files:
        if file.endswith('_genomic.fna') and not file.endswith('from_genomic.fna'):
            to_copy = file
            break

    if to_copy:
        # copy the file to the output dir
        shutil.copy(to_copy, os.path.join(output_dir, f"reference.fna"))
        full_reference_name = os.path.basename(to_copy)
    else:
        raise FileNotFoundError(f"Could not find the reference genome for {reference}")

    # do the same thing for the simulation
    simulation_files = glob.glob(os.path.join(genome_dir, "**", f"{simulation}*"), recursive=True)
    to_copy = ''
    for file in simulation_files:
        if file.endswith('_genomic.fna') and not file.endswith('from_genomic.fna'):
            to_copy = file
            break

    if to_copy:
        # copy the file to the output dir
        shutil.copy(to_copy, os.path.join(output_dir, f"simulation.fna"))
        shutil.copy(to_copy, os.path.join(output_dir, f"simulated_mg.fq"))
        full_simulation_name = os.path.basename(to_copy)
    else:
        raise FileNotFoundError(f"Could not find the reference genome for {simulation}")

    # write the distance to a file
    with open(os.path.join(output_dir, "distance.txt"), 'w') as f:
        f.write(f"ANI: {distance}")

    # write the fake simulation_counts.csv file
    with open(os.path.join(output_dir, "simulation_counts.csv"), 'w') as f:
        f.write(f"{full_simulation_name},1")

    # make the picklist
    with open(os.path.join(output_dir, "known_names_picklist.txt"), 'w') as f:
        f.write(f"internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n")
    found_it = False
    with open("MANIFEST.csv", 'r') as f:
        lines = f.readlines()
        for line in lines:
            if len(line.split(',')) >= 2:
                name = line.split(',')[-2]
                if name.startswith(reference):
                    with open(os.path.join(output_dir, "known_names_picklist.txt"), 'a') as f:
                        f.write(line)
                    found_it = True
                    break
    if not found_it:
        raise FileNotFoundError(f"Could not find the reference genome for {reference} in the manifest. DO NOT USE "
                                f"THESE RESULTS!")


if __name__ == '__main__':
    main()

