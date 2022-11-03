#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Make the sourmash sketch
sourmash sketch dna -f -p k=31,scaled=1000,abund -o simulated_mg.fq.sig simulated_mg.fq
# Run gather
sourmash gather --dna --threshold-bp 100 simulated_mg.fq.sig without_unknown_db.sig -o gather_results.csv

# then run our approach
python ../../ref_matrix.py --ref_file without_unknown_db.sig  --ksize 31 --out_prefix default_EU
