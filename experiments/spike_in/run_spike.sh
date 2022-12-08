#!/usr/bin/env bash
set -e
set -u
set -o pipefail

# This script will prepare for the spike-in experiment

# Train our method on everything
cd ..
python ../../ref_matrix.py --ref_file ../formatted_db.sig --ksize 31 --out_prefix formatted_db_ --max_thresh 5
cd spike_in

# Sketch the real metagenome
sourmash sketch dna -f -p k=31,scaled=1000,abund -o 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq

# Run sourmash gather
sourmash gather --dna --threshold-bp 0 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig ../formatted_db.sig -o gather_results.csv

# get the names of the detected organisms
cut -d',' -f 10 gather_results.csv | grep -v name > detected_orgs.txt

# get the names of the undetected organisms
grep -v -f detected_orgs.txt ../MANIFEST.csv | cut -d',' -f 10 > absent_names.txt

#iter=1
coverageValues=(".5" ".25" ".125" ".0625" ".03125" ".015625" ".0078125" ".00390625" ".001953125" ".0009765625") 
for iter in `seq 1 100`; do
# For commenting out stuff I've already done :<<'END' END
# Randomly pull one of these out
cat absent_names.txt | shuf | head -n1 > spike_in_name.txt

# Find the reference and copy it here
spikeName=$(find .. -name `cat spike_in_name.txt`)
cp ${spikeName} to_spike_in.fna

# Subsample to various coverage amounts
# This assumes that the repo KEGG_sketching_annotation is adjacent to Estimating_Unknowns
for cov in ${coverageValues[@]}
do
	../../../KEGG_sketching_annotation/utils/bbmap/./randomreads.sh ref=to_spike_in.fna overwrite=t out=coverage_${cov}X.fna coverage=0${cov}
done

# Sketch all of these
for cov in ${coverageValues[@]}
do 
	sourmash sketch dna -f -p k=31,scaled=1000,abund -o coverage_${cov}X.fna.sig coverage_${cov}X.fna
done

# Merge these into the sketch of the real metagenome
for cov in ${coverageValues[@]}
do 
	sourmash sig merge -o merged_${cov}X.sig coverage_${cov}X.fna.sig 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig
done

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

done
