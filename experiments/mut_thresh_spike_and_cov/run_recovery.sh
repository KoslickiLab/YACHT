#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Run our abundance recovery
for sampleCov in ${coverageValues[@]}
do
        for covThresh in ${coverageValues[@]}
        do
                echo "python ../../recover_abundance.py --ref_file ../formatted_db_ref_matrix_processed.npz --ksize 31 --sample_file merged_${sampleCov}X.sig --outfile EU_results_sample_cov_${sampleCov}X_cov_thresh_${covThresh}X.csv --min_coverage ${covThresh}"
        done
done | parallel -j 50

#BACK_PID=$!
#wait $BACK_PID

# Prep the output CSV file
echo "sample_coverage, coverage_threshold, recovered_kmer_abundance" > results_iter_${iter}.csv

# Get the results
for sampleCov in ${coverageValues[@]}
do
        for covThresh in ${coverageValues[@]}
        do
                relAb=$(grep -f spike_in_name.txt EU_results_sample_cov_${sampleCov}X_cov_thresh_${covThresh}X.csv | cut -d',' -f17)
                echo "${sampleCov}, ${covThresh}, ${relAb}" >> results_iter_${iter}.csv
        done
done

