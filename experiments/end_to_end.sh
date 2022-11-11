#!/usr/bin/env bash
set -e
set -u
set -o pipefail
# This first command need only be run once
echo "Making database"
#python ../../KEGG_sketching_annotation/scripts/create_genome_ref_db.py ./ref_genomes_3/reference_genomes 1046_database 1046

numReads=10000000

simsFolder=sims-$(date +"%s")
mkdir ${simsFolder}
echo "Creating simulation"
cd ${simsFolder}
ln -s ../formatted_db.fasta formatted_db.fasta
cd ..
python run_sim.py --genomes_folder ref_genomes_3/reference_genomes --out_folder ${simsFolder} --num_reads ${numReads} --num_orgs 1000


# Note: I could do a loop here if I felt like it
N=500
# remove N genomes from the training datbase
echo "Making ${N} genomes unkown"
cat ${simsFolder}/simulation_counts.csv | shuf | head -n ${N} | cut -d',' -f1 | sort > ${simsFolder}/unknown_names.txt
cat ${simsFolder}/simulation_counts.csv | cut -d',' -f1 | sort > ${simsFolder}/all_names.txt
comm -2 -3 ${simsFolder}/all_names.txt ${simsFolder}/unknown_names.txt > ${simsFolder}/known_names.txt
sed -n 2p MANIFEST.csv > ${simsFolder}/known_names_picklist.txt
cat ${simsFolder}/known_names.txt |cut -d'.' -f1 | xargs -I{} grep {} MANIFEST.csv >> ${simsFolder}/known_names_picklist.txt
#echo "name" > ${simsFolder}/known_names_picklist.txt
#cat ${simsFolder}/known_names.txt >> ${simsFolder}/known_names_picklist.txt

# remove those from the training datbase
echo "Removing them from the ref db"
#../../KEGG_sketching_annotation/utils/bbmap/./filterbyname.sh in=${simsFolder}/formatted_db.fasta out=${simsFolder}/without_unknown_db.fasta names=${simsFolder}/unknown_names.txt include=f overwrite=true
# Sketching the reference 
echo "Sketching the reference"
#sourmash sketch dna -f -p k=31,scaled=1000,abund -o ${simsFolder}/without_unknown_db.sig --singleton ${simsFolder}/without_unknown_db.fasta
sourmash sig extract --picklist ${simsFolder}/known_names_picklist.txt:md5:md5 formatted_db.sig -o ${simsFolder}/without_unknown_db.sig 
echo "Making the EU dictionary"
python ../ref_matrix.py --ref_file ${simsFolder}/without_unknown_db.sig  --ksize 31 --out_prefix ${simsFolder}/default_EU_
# then run the methods
echo "sketching the simulation"
sourmash sketch dna -f -p k=31,scaled=1000,abund -o ${simsFolder}/simulated_mg.fq.sig ${simsFolder}/simulated_mg.fq
# Run gather
echo "running gather"
sourmash gather --dna --threshold-bp 100 ${simsFolder}/simulated_mg.fq.sig ${simsFolder}/without_unknown_db.sig -o ${simsFolder}/gather_results.csv

# then run our approach
python ../recover_abundance.py --ref_file ${simsFolder}/default_EU_ref_matrix_processed.npz  --ksize 31 --sample_file ${simsFolder}/simulated_mg.fq.sig --outfile ${simsFolder}/EU_results_default.csv

numUnknownReads=$(grep -v -f ${simsFolder}/known_names.txt ${simsFolder}/simulation_counts.csv | cut -d',' -f2 | awk '{SUM+=$1}END{print SUM}')
unknownPercent=$(echo "scale=8; ${numUnknownReads} / ${numReads}" | bc)
sourmashEstimateKnown=$(cut -d',' -f 5 ${simsFolder}/gather_results.csv | tail -n +2 | awk '{SUM+=$1}END{print SUM}')
sourmashEstimateUnknown=$(echo "scale=8; 1 - ${sourmashEstimateKnown}" | bc)
ourEstimate=$(sed -n 2p ${simsFolder}/EU_results_default.csv | rev | cut -d',' -f1 | rev)
echo "true_unknown,sourmash_estimate,our_estimate" > ${simsFolder}/results.txt
echo "${unknownPercent},${sourmashEstimateUnknown},${ourEstimate}" >> ${simsFolder}/results.txt
echo "${unknownPercent},${sourmashEstimateUnknown},${ourEstimate}" >> all_results.txt
rm ${simsFolder}/simulated_mg.fq
