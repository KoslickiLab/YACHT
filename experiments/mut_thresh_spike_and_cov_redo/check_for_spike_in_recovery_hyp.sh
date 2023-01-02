#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail

# For each of the experiments, find who they are supposed to be similar to and check if it's in the sample
# prep the output csv
echo "start"
echo "spike_coverage@coverage_threshold@max_ani@rel_ab@mut_thresh@spike_md5short@spike_name@match_name@num_exclusive_kmers_with_coverage@num_matches@acceptance_threshold_with_coverage@FP_flag@FN_flag" > results.csv
coverageValues=("1" "0.1" "0.01" "0.001")
mutThreshes=("0.05" "0.01" "0.001")
for mutThresh in "${mutThreshes[@]}"
do
	for spikeCov in "${coverageValues[@]}"
	do
        	for covThresh in "${coverageValues[@]}"
        	do
        	        numLines=$(wc -l sigs_md5_to_accession_to_gtdb_location_mut_${mutThresh}.txt | cut -d' ' -f1)
        	        it=0
                	for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${mutThresh}.txt`
                	do
                	        #increment the iterator
                	        it=$((it+1))
                	        # if the iterator is a multiple of 100, print a message
                	        if [ $((it%100)) -eq 0 ]; then
                	                echo "On line ${it} of ${numLines}"
                	        fi
                        	md5short=$(echo ${line} | cut -d',' -f2)
                          md5long=$(echo ${line} | cut -d',' -f1)
                          match=$(grep -P -m 1 ${md5long} in_gtdb_similar_to_EU_not_in_sample_mut_${mutThresh}_@.tsv)
                          # Just take the accession, so we don't need to deal with commas in names
                          gtdbName=$(echo $match | cut -d'@' -f1 | cut -d' ' -f1)
                          gtdbMD5=$(echo $match | cut -d'@' -f2)
                          maxANI=$(echo $match | cut -d'@' -f3)
                          matchingOrg=$(echo $match | cut -d'@' -f4)
                          gtdbShortName=$(echo ${gtdbName} | cut -d'.' -f1)
                          line=$(grep ${matchingOrg} EU_on_spikes_cov_${spikeCov}_hyp/${md5short}_${spikeCov}X_cov_thresh_${covThresh}X_mut_thresh_${mutThresh}_hyp.csv)
                          if [ $? -eq 0 ]; then
                            relAb=$(echo ${line} | cut -d',' -f 13)
                            num_exclusive_kmers_with_coverage=$(echo ${line} | cut -d',' -f 15)
                            num_matches=$(echo ${line} | cut -d',' -f 16)
                            acceptance_threshold_with_coverage=$(echo ${line} | cut -d',' -f 18)
                            # if relAb is 0, but maxANI is greater than mutThresh, then it's a false negative
                            if [ "$(echo "${relAb} == 0.0" | bc)" -eq 1 ] && [ "$(echo "1-${maxANI} < ${mutThresh}" | bc)" -eq 1 ]; then
                              FN_flag=1
                            else
                              FN_flag=0
                            fi
                            # if relAb is greater than 0, but maxANI is less than mutThresh, then it's a false positive
                            if [ "$(echo "${relAb} > 0.0" | bc)" -eq 1 ] && [ "$(echo "1-${maxANI} > ${mutThresh}" | bc)" -eq 1 ]; then
                              FP_flag=1
                            else
                              FP_flag=0
                            fi
                          else
                            relAb="NaN"
                            num_exclusive_kmers_with_coverage="NaN"
                            num_matches="NaN"
                            acceptance_threshold_with_coverage="NaN"
                          fi

                          echo "${spikeCov}@${covThresh}@${maxANI}@${relAb}@${mutThresh}@${md5short}@${gtdbName}@${matchingOrg}@${num_exclusive_kmers_with_coverage}@${num_matches}@${acceptance_threshold_with_coverage}@${FP_flag}@${FN_flag}" >> results.csv
			done
		done
	 done
done
echo "finish"
