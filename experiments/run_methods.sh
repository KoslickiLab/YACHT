#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Make the sourmash sketch
#sourmash sketch dna -f -p k=31,scaled=1000,abund -o simulated_mg_10M.fq.sig simulated_mg_10M.fq
# Run gather
sourmash gather --dna --threshold-bp 100 simulated_mg_10M.fq.sig without_unknown_db.sig -o gather_results_10M.csv

# then run our approach
python ../../recover_abundance.py --ref_file default_EUref_matrix_processed.npz  --ksize 31 --sample_file simulated_mg_10M.fq.sig --outfile EU_results_default_10M.csv
