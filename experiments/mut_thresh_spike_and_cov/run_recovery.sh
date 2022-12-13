#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Run our abundance recovery
mutThreshs=("0.001" "0.005" "0.01" "0.02" "0.05" "0.1" "0.15")
coverageValues=(".25" ".0625" ".015625" ".00390625" ".0009765625")
for mut in ${mutThreshs[@]}
do
	for spikeCov in ${coverageValues[@]}
	do
        	for covThresh in ${coverageValues[@]}
		do
			for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location.txt`
	        	do
               			md5short=$(echo ${line} | cut -d',' -f2)
                		echo "python ../../recover_abundance.py --ref_file ../formatted_db_ref_matrix_processed.npz --ksize 31 --sample_file spikes_cov_${spikeCov}/${md5short}_spiked_sample.sig.zip --outfile EU_on_spikes_cov_${spikeCov}/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X_mut_thresh_${mut}.csv --min_coverage ${covThresh}" --mut_thresh ${mut}
	        	done | parallel -j 50
		done
	done
done

