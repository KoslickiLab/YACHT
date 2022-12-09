#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# For each of the experiments, find who they are supposed to be similar to and check if it's in the sample
# prep the output csv
echo "spike_coverage,coverage_threshold,max_ani,rel_ab" > results.csv
coverageValues=(".25" ".0625" ".015625" ".00390625" ".0009765625")
for spikeCov in ${coverageValues[@]}
do
        for covThresh in ${coverageValues[@]}
        do
                for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location.txt`
                do
                        md5short=$(echo ${line} | cut -d',' -f2)
			md5long=$(echo ${line} | cut -d',' -f1)
			read -r gtdbName gtdbMD5 maxANI matchingOrg <<< $(grep -m 1 ${md5long} in_gtdb_similar_to_EU_not_in_sample.tsv)
			gtdbShortName=$(echo ${gtdbName} | cut -d'.' -f1)
			relAb=$(grep ${gtdbShortName} EU_on_spikes_cov_${spikeCov}/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X.csv)
			echo "${spikeCov},${covThresh},${maxANI},${relAb}" >> results.csv
		done
	done
done

