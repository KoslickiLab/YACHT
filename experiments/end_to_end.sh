#!/usr/bin/env bash
set -e
set -u
set -o pipefail
# This first command need only be run once
echo "Making database"
# python ../../KEGG_sketching_annotation/scripts/create_genome_ref_db.py ./ref_genomes_3/reference_genomes 1046_database 1046

numReads=10000000

simsFolder=sims-$(date +"%s")
mkdir ${simsFolder}
echo "Creating simulation"
cd ${simsFolder}
ln -s ../formatted_db.fasta formatted_db.fasta
python ../run_sim.py --genomes_folder ../ref_genomes_3/reference_genomes --out_folder . --num_reads ${numReads} --num_orgs 200
echo "sketching the simulation"
sourmash sketch dna -f -p k=31,scaled=1000,abund -o simulated_mg.fq.sig simulated_mg.fq
rm simulated_mg.fq
rm -rf ref

# Note: I could do a loop here if I felt like it
N=100
# remove N genomes from the training datbase
echo "Making ${N} genomes unkown"
cat simulation_counts.csv | shuf | head -n ${N} | cut -d',' -f1 | sort > unknown_names.txt
cat simulation_counts.csv | cut -d',' -f1 | sort > all_names.txt
# This is if you want true negatives
grep '>' formatted_db.fasta | sed 's/>//g' | sort > all_names.txt
comm -2 -3 all_names.txt unknown_names.txt > known_names.txt
sed -n 2p ../MANIFEST.csv > known_names_picklist.txt
cat known_names.txt | cut -d'.' -f1 | xargs -I{} grep {} ../MANIFEST.csv >> known_names_picklist.txt

# remove those from the training datbase
echo "Removing them from the ref db"
# Sketching the reference 
echo "Sketching the reference"
sourmash sig extract --picklist known_names_picklist.txt:md5:md5 ../formatted_db.sig -o without_unknown_db.sig 
echo "Making the EU dictionary"
python ../../ref_matrix.py --ref_file without_unknown_db.sig  --ksize 31 --out_prefix default_EU_ --max_thresh 5
# then run the methods
# Run gather
echo "running gather"
sourmash gather --dna --threshold-bp 100 simulated_mg.fq.sig without_unknown_db.sig -o gather_results.csv #--output-unassigned gather_unassigned.sig

# then run our approach
python ../../recover_abundance.py --ref_file default_EU_ref_matrix_processed.npz  --ksize 31 --sample_file simulated_mg.fq.sig --outfile EU_results_default.csv --p_val 0.01

#numUnknownReads=$(grep -v -f known_names.txt simulation_counts.csv | cut -d',' -f2 | awk '{SUM+=$1}END{print SUM}')
#numUnknownReads=$(cat unknown_names.txt | xargs -I{} grep {} simulation_counts.csv | cut -d',' -f2 | awk '{SUM+=$1}END{print SUM}')
numUnknownReads=$(../calculate_unknown_percent.py -d . -r ../formatted_db.sig)
unknownPercent=$(echo "scale=8; ${numUnknownReads} / ${numReads}" | bc)
sourmashEstimateKnown=$(cut -d',' -f 5 gather_results.csv | tail -n +2 | awk '{SUM+=$1}END{print SUM}')
sourmashEstimateUnknown=$(echo "scale=8; 1 - ${sourmashEstimateKnown}" | bc)
ourEstimate=$(sed -n 2p EU_results_default.csv | rev | cut -d',' -f1 | rev)
echo "true_unknown,sourmash_estimate,our_estimate,num_reads" > results.txt
echo "${unknownPercent},${sourmashEstimateUnknown},${ourEstimate},${numReads}" >> results.txt
echo "${unknownPercent},${sourmashEstimateUnknown},${ourEstimate},${numReads}" >> ../all_results.txt
