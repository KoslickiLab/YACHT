#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# This assumes that you have run the run_sim.py file and put it in a directory `sims`. It also assumes that the KEGG_sketching repo is next to the Estimating_Unknowns repo

N=10
# remove N genomes from the training datbase
cat sims/simulation_counts.csv | shuf | head -n ${N} | cut -d',' -f1 > sims/unknown_names.txt

# remove those from the training datbase
../../../KEGG_sketching_annotation/utils/bbmap/./filterbyname.sh in=formatted_db.fasta out=without_unknown_db.fasta names=unknown_names.txt include=f
