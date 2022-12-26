#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# Run our abundance recovery
#mutThreshs=("0.1" "0.05" "0.01" "0.001")
mutThreshs=("0.05" "0.01" "0.001")
minMut=0.001
coverageValues=("1" "0.1" "0.01" "0.001")
for covThresh in "${coverageValues[@]}"
do
	for mut in "${mutThreshs[@]}"
	do
        	for spikeCov in "${coverageValues[@]}"
		do
			mkdir -p EU_on_spikes_cov_${spikeCov}_hyp
			for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${mut}.txt`
	        	do
               			md5short=$(echo ${line} | cut -d',' -f2)
                		echo "python ../../recover_abundance.py --ref_file formatted_db_mut_${mut}_ref_matrix_processed.npz --ksize 31 --sample_file spikes_cov_${spikeCov}/${md5short}_spiked_sample.sig.zip --outfile EU_on_spikes_cov_${spikeCov}_hyp/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X_mut_thresh_${mut}_hyp.csv --min_coverage ${covThresh}" --mut_thresh ${mut} --recovery_method h
	        	done | parallel -j 200
		done
	done
done


