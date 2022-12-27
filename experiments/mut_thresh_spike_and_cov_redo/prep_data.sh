#!/usr/bin/env bash
set -e
set -u
set -o pipefail
#: <<'END'

# These are the mutation rates I'll be using:
mutThreshs=("0.001" "0.01" "0.05" "0.1")
minMut=0.001
coverageValues=("1" "0.1" "0.01" "0.001")
# and these are the coverage values
: <<'END'
# This script will get the required data for the spike-in experiment

# sketch the reference database
sourmash sketch dna -p k=31,scaled=1000,abund -o formatted_db.sig --singleton ../formatted_db.fasta

# Get the GTDB sourmash database
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip

# extract the manifest
sourmash sig manifest gtdb-rs207.genomic-reps.dna.k31.zip -o MANIFEST.csv --no-rebuild-manifest
#END
# Make the required folders
mkdir -p sigs
mkdir -p spikes
mkdir -p EU_on_spikes


# get the real metagenome
# FIXME: will want to make this reach directly out to NCBI
# cp /data/shared_data/TwinsStudy/data/MZ/36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.gz .
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR911/ERR911999/36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.gz

# sketch the metagenome
sourmash sketch dna -p k=31,scaled=1000,abund 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.gz -o 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig

# Run sourmash gather
sourmash gather --dna --threshold-bp 0 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig formatted_db.sig -o gather_results.csv

# get the names of the detected organisms
cut -d',' -f 10 gather_results.csv | grep -v name > detected_orgs.txt
END
# get the names of the undetected organisms
sourmash sig collect formatted_db.sig -o MANIFEST.csv -F csv --merge-previous
grep -v -f detected_orgs.txt ../MANIFEST.csv | cut -d',' -f 10 > absent_names.txt

# find organisms in GTDB that are similar to the reference database, but not in the sample
# find organisms in GTDB that are in the sample
sourmash gather --dna --threshold-bp 0 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig gtdb-rs207.genomic-reps.dna.k31.zip -o gather_on_gtdb.csv --save-prefetch-csv gather_on_gtdb_prefetch.csv --save-prefetch gather_on_gtdb_prefetch.sig
# from the MANIFEST.csv, remove those that are in gather_on_gtdb.csv
comm -3 -2 <(cat MANIFEST.csv | cut -d',' -f2 | sort) <(cat gather_on_gtdb.csv | cut -d',' -f12 | sort) | grep -v "#" > md5_in_gtdb_not_in_sample.txt
# add the md5 header
cat <(echo md5) <(cat md5_in_gtdb_not_in_sample.txt) > md5_in_gtdb_not_in_sample.csv
# extract these from the GTDB database
sourmash sig extract --picklist md5_in_gtdb_not_in_sample.csv:md5:md5 -o gtdb-rs207.genomic-reps.dna.k31_not_in_sample.sig gtdb-rs207.genomic-reps.dna.k31.zip

# sketch the whole reference database in order to prep for a prefilter search
sourmash sketch dna -p k=31,scaled=1000,abund ../formatted_db.fasta -o formatted_db_merged.sig

#python ../../ref_matrix.py --ref_file formatted_db.sig --ksize 31 --out_prefix formatted_db_mut_${mut} --max_thresh 5 --mut_thresh 0.05
# FIXME: will need to use the merged db, but only consisting of those that survive the ref_matrix build process
#sourmash sig merge -o formatted_db_merged.sig formatted_db

# sourmash prefetch seems really slow for som reason, so let's use a gather and save the prefetch stuff, no need to worry about --threshold-bp since we won't actually use the gather results here, just the prefetch
# This will take a long time
sourmash gather --dna formatted_db_merged.sig gtdb-rs207.genomic-reps.dna.k31_not_in_sample.sig -o gather_formatted_db_merged_on_gtdb_not_in_sample.csv --save-prefetch-csv gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.csv --save-prefetch gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig

# the gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig file contains candidates from GTDB that are similar to the reference, but not in the sample. So now we need to:
# extract these candidate ones
sourmash sig extract --picklist gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.csv:match_md5:md5short -o gtdb-rs207.genomic-reps.dna.k31_not_in_sample_similar_to_merged_reference.sig gtdb-rs207.genomic-reps.dna.k31.zip

# build all of the references
# loop through the mutation rates and make the dictionary
for mut in "${mutThreshs[@]}"
do
  python ../../ref_matrix.py --ref_file formatted_db.sig --ksize 31 --out_prefix formatted_db_mut_${mut}_ --mut_thresh ${mut}
done

# then calculate the actual max ANI to the reference
for mut in "${mutThreshs[@]}"
do
  python calculate_ANI_of_gtdb_to_reference.py --mutation_rate ${mut} --gtdb gtdb-rs207.genomic-reps.dna.k31_not_in_sample_similar_to_merged_reference.sig --reference_database_full formatted_db.sig
done
#END
for cov in "${coverageValues[@]}"
do
        mkdir -p sigs_cov_${cov}
        mkdir -p sigs_cov_${cov}/reads
        mkdir -p spikes_cov_${cov}
        mkdir -p EU_on_spikes_cov_${cov}
done

# Will need to do this for each mutation value, since I'll need the sigs_md5_to_accession_to_gtdb_location.txt files to properly consider only those that are actually similar to something in the dictionary
for mut in "${mutThreshs[@]}"
do
    echo "Start loop"
    # then extract all the relevant sigs
    cut -f2 in_gtdb_similar_to_EU_not_in_sample_mut_${mut}.tsv > in_gtdb_similar_to_EU_not_in_sample_mut_${mut}_md5.txt

    sourmash sig split -f --output-dir sigs_mut_${mut} gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig --picklist in_gtdb_similar_to_EU_not_in_sample_mut_${mut}_md5.txt:gtdb_md5:md5

    # Get the accessions of each of the spikes
    sourmash sig collect -F csv --merge-previous -o sigs_manifest_mut_${mut}.csv sigs_mut_${mut}/*
    cut -d',' -f2,3,10 sigs_manifest_mut_${mut}.csv | sed 's/"//g' | cut -d' ' -f1 | grep -v \# > sigs_md5_to_accession_mut_${mut}.txt

    # find the location of the real GTDB genomes
    #cut -d',' -f3 sigs_md5_to_accession_mut_${mut}.txt | tail -n +2 | xargs -P 10 -I{} grep -m1 {} /data/shared_data/GTDB/gtdb_genomes_reps_r207/file_list.txt > sigs_gtdb_file_locations_mut_${mut}.txt
    cut -d',' -f3 sigs_md5_to_accession_mut_${mut}.txt | tail -n +2 | parallel -j 10 --keep-order grep -m1 {} /data/shared_data/GTDB/gtdb_genomes_reps_r207/file_list.txt > sigs_gtdb_file_locations_mut_${mut}.txt

    # add these to the accession list
    paste -d',' <(cat sigs_md5_to_accession_mut_${mut}.txt) <(sed '1s/^/gtdb_file_location\n/' sigs_gtdb_file_locations_mut_${mut}.txt) > sigs_md5_to_accession_to_gtdb_location_mut_${mut}.txt
    echo "End loop"
done

# change the delimiter to @
for mut in "${mutThreshs[@]}"
do
    sed 's/\t/@/g' in_gtdb_similar_to_EU_not_in_sample_mut_${mut}.tsv > in_gtdb_similar_to_EU_not_in_sample_mut_${mut}_@.tsv
done





# For each of these GTDB genomes, get X coverage of the genome, sketch it, and then stick it in a folder

# do this for the min mutation rate, as this will contain the most genomes
# reduce the coverage
for cov in "${coverageValues[@]}"
do
        for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${minMut}.txt`
        do
                md5short=$(echo ${line} | cut -d',' -f2)
                fileLoc=$(echo ${line} | cut -d',' -f4)
                # bbmap to coverage amount
                #../../../KEGG_sketching_annotation/utils/bbmap/./randomreads.sh ref=${fileLoc} overwrite=t out=sigs_cov_${cov}/reads/${md5short}.fna coverage=0${cov}
                tempDir=$(mktemp -d)
                echo "cd $tempDir; /data/dmk333/Repositories/KEGG_sketching_annotation/utils/bbmap/./randomreads.sh  ref=${fileLoc} overwrite=t out=/data/dmk333/Repositories/Estimating_Unknowns/experiments/mut_thresh_spike_and_cov_redo/sigs_cov_${cov}/reads/${md5short}.fna coverage=0${cov}"
        done | parallel -j 200
done

# then sketch each one
for cov in "${coverageValues[@]}"
do
        for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${minMut}.txt`
        do
                md5short=$(echo ${line} | cut -d',' -f2)
                fileLoc=$(echo ${line} | cut -d',' -f4)
                echo sourmash sketch dna -f -p k=31,scaled=1000,abund -o sigs_cov_${cov}/${md5short}.k=31.scaled=1000.DNA.dup=0.63.sig.zip sigs_cov_${cov}/reads/${md5short}.fna
        done | parallel -j 50
done

# then create all the spiked samples
for cov in "${coverageValues[@]}"
do
        for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location_mut_${minMut}.txt`
        do
                md5short=$(echo ${line} | cut -d',' -f2)
                fileLoc=$(echo ${line} | cut -d',' -f4)
                echo sourmash sig merge -o spikes_cov_${cov}/${md5short}_spiked_sample.sig.zip sigs_cov_${cov}/${md5short}.k=31.scaled=1000.DNA.dup=0.63.sig.zip 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig
        done | parallel -j 50
done
