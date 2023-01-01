#!/usr/bin/env bash
#set -e
#set -u
#set -o pipefail
# This first command need only be run once
echo "Making database"
# python ../../KEGG_sketching_annotation/scripts/create_genome_ref_db.py ./ref_genomes_3/reference_genomes 1046_database 1046

coverageValues=("1" "0.1" "0.01" "0.001")
#coverage=.1
for coverage in "${coverageValues[@]}"
do
numReads=10000000

simsFolder=sims-$(date +"%s")-coverage-${coverage}
mkdir ${simsFolder}
echo "Creating simulation"
cd ${simsFolder}
ln -s ../formatted_db.fasta formatted_db.fasta
python ../run_sim.py --genomes_folder ../ref_genomes_3/reference_genomes --out_folder . --num_reads ${numReads} --num_orgs 200 #--uniform
echo "sketching the simulation"
sourmash sketch dna -f -p k=31,scaled=1000,abund -o simulated_mg.fq.sig simulated_mg.fq
rm simulated_mg.fq
rm -rf ref

# Note: I could do a loop here if I felt like it
N=100
# remove N genomes from the training datbase
echo "Making ${N} genomes unkown"
cat simulation_counts.csv | shuf | head -n ${N} | cut -d',' -f1 | sort > unknown_names.txt  # randomly select N names that showed up in the simulation and mark them as unknown
cat simulation_counts.csv | cut -d',' -f1 | sort > all_names_in_sim.txt  # get all the orgs that actually showed up in the simulation
grep '>' formatted_db.fasta | sed 's/>//g' | sort > all_names.txt  # get all the names of orgs we downloaded
comm -2 -3 all_names.txt unknown_names.txt > known_names.txt  # subtract from all the names the ones marked as unknown. This will form the reference database for both methods
sed -n 2p ../MANIFEST.csv > known_names_picklist.txt  # get the lines from the manifest in order to pick them out of the pre-computed sketches
cat known_names.txt | cut -d'.' -f1 | xargs -I{} grep {} ../MANIFEST.csv >> known_names_picklist.txt

# What I should really be doing is taking the WHOLE reference database, removing the unknown, and then using that

# remove those from the training datbase
echo "Removing them from the ref db"
# Sketching the reference 
echo "Sketching the reference"
sourmash sig extract --picklist known_names_picklist.txt:md5:md5 ../formatted_db.sig -o without_unknown_db.sig # all organisms (in the simulation or not) that are not marked as unknown
echo "Making the EU dictionary"
python ../../ref_matrix.py --ref_file without_unknown_db.sig --ksize 31 --out_prefix default_EU_
# then run the methods
# Run gather
echo "running gather"
sourmash gather --dna --threshold-bp 100 simulated_mg.fq.sig without_unknown_db.sig -o gather_results.csv #--output-unassigned gather_unassigned.sig

# then run our approach
python ../../recover_abundance.py --ref_file default_EU_ref_matrix_processed.npz  --ksize 31 --sample_file simulated_mg.fq.sig --outfile EU_results_default.csv --min_coverage ${coverage} --recovery_method h

# also print the binary stats
../calculate_unknown_percent.py -d . -r ../formatted_db.sig
python ../calculate_binary_stats_hyp.py -d .
cd ..
done
