#!/bin/bash
set -e
set -u
set -o pipefail

# for a bunch of mut_thresh values, run our method
mutThreshs=("0.001" "0.005" "0.01" "0.02" "0.05" "0.1" "0.15" "0.2" "0.3")
for mut in ${mutThreshs[@]}
do
	for spike in `ls spikes/*`
	do
		name=$(basename ${spike})
		echo python ../../recover_abundance.py --ref_file ../formatted_db_ref_matrix_processed.npz --ksize 31 --sample_file ${spike} --outfile EU_on_spikes/EU_results_mut_thresh_${mut}_${name}.csv --mut_thresh ${mut}
	done
done | parallel -j 200
