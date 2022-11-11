#!/usr/bin/env python
# This script will take all the organisms marked as unknown and check if there is something in the known
# part that is 0.95 ANI or higher. If so, it will mark it as actually known.
import os
import argparse
import sourmash


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Adds edge lengths to the KEGG hierarchy')
    parser.add_argument('-d', '--dir', help='Directory containing the unknown and known signatures', required=True)
    parser.add_argument('-r', '--ref', help='Reference/training database sketches/signature', required=True)
    args = parser.parse_args()
    sims_dir = args.dir
    ref_db = args.ref


    if not os.path.exists(sims_dir):
        raise FileNotFoundError(f"Could not find {sims_dir}")
    if not os.path.exists(ref_db):
        raise FileNotFoundError(f"Could not find {ref_db} in the same directory as where this was run")
    all_names_file = os.path.join(sims_dir, 'all_names.txt')
    if not os.path.exists(all_names_file):
        raise FileNotFoundError(f"Could not find {all_names_file} in the same directory as {sims_dir}")
    unknown_names_file = os.path.join(sims_dir, 'unknown_names.txt')
    if not os.path.exists(unknown_names_file):
        raise FileNotFoundError(f"Could not find {unknown_names_file} in the same directory as {sims_dir}")
    known_names_file = os.path.join(sims_dir, 'known_names.txt')
    if not os.path.exists(known_names_file):
        raise FileNotFoundError(f"Could not find {known_names_file} in the same directory as {sims_dir}")
    simulation_counts_file = os.path.join(sims_dir, 'simulation_counts.csv')
    if not os.path.exists(simulation_counts_file):
        raise FileNotFoundError(f"Could not find {simulation_counts_file} in the same directory as {sims_dir}")

    sketches = list(sourmash.load_file_as_signatures(ref_db))
    # Get all the unknown names
    with open(unknown_names_file, 'r') as f:
        unknown_names = [line.strip().split('.')[0] for line in f.readlines() if line.strip()]
    # Get all the known names
    with open(known_names_file, 'r') as f:
        known_names = [line.strip().split('.')[0] for line in f.readlines() if line.strip()]

    name_2_sketch = {sketch.name.split('.')[0]: sketch for sketch in sketches}
    # iterate through the unknown names, check if there is a known name that is 0.95 or higher
    actually_unknown = []
    for unknown_name in unknown_names:
        is_similar_to_known = False
        unknown_sketch = name_2_sketch[unknown_name]
        for known_name in known_names:
            known_sketch = name_2_sketch[known_name]
            ani = unknown_sketch.max_containment_ani(known_sketch, estimate_ci=True).ani_high
            #ani = unknown_sketch.max_containment_ani(known_sketch).ani
            if ani and ani >= 0.95:
                is_similar_to_known = True
                print(f"{unknown_name} is similar to {known_name} with ANI {ani}")
                break
        if not is_similar_to_known:
            actually_unknown.append(unknown_name)

    # now sum up the counts for the actually unknowns
    name_2_counts = {}
    with open(simulation_counts_file, 'r') as f:
        for line in f.readlines():
            name, count = line.strip().split(',')
            name = name.split('.')[0]
            name_2_counts[name] = int(count)
    unknown_counts = [name_2_counts[name] for name in actually_unknown]
    print(f"{sum(unknown_counts)}")


if __name__ == '__main__':
    main()
