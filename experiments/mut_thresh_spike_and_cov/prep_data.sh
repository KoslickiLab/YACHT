#!/usr/bin/env bash
set -e
set -u
set -o pipefail
: <<'END'
# This script will get the required data for the spike-in experiment

# Get the GTDB sourmash database
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip

# extract the manifest
sourmash sig manifest gtdb-rs207.genomic-reps.dna.k31.zip -o MANIFEST.csv --no-rebuild-manifest

# Make the required folders
mkdir sigs
mkdir spikes
mkdir EU_on_spikes

# get the real metagenome
cp /data/shared_data/TwinsStudy/data/MZ/36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.gz .

# sketch the metagenome
sourmash sketch dna -p k=31,scaled=1000,abund 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.gz -o 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig

# Run sourmash gather
sourmash gather --dna --threshold-bp 0 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig ../formatted_db.sig -o gather_results.csv

# get the names of the detected organisms
cut -d',' -f 10 gather_results.csv | grep -v name > detected_orgs.txt

# get the names of the undetected organisms
grep -v -f detected_orgs.txt ../MANIFEST.csv | cut -d',' -f 10 > absent_names.txt

# find organisms in GTDB that are similar to the reference database, but not in the sample
# find organisms in GTDB that are in the sample
sourmash gather --dna --threshold-bp 0 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig gtdb-rs207.genomic-reps.dna.k31.zip -o gather_on_gtdb.csv --save-prefetch-csv gather_on_gtdb_prefetch.csv --save-prefetch gather_on_gtdb_prefetch.sig
# from the MANIFEST.csv, remove those that are in gather_on_gtdb.csv
comm -3 -2 <(cat MANIFEST.csv | cut -d',' -f2 | sort) <(cat gather_on_gtdb.csv | cut -d',' -f12 | sort) | grep -v "#" > md5_in_gtdb_not_in_sample.txt# add the md5 header
cat <(echo md5) <(cat md5_in_gtdb_not_in_sample.txt) > md5_in_gtdb_not_in_sample.csv
# extract these from the GTDB database
sourmash sig extract --picklist md5_in_gtdb_not_in_sample.csv:md5:md5 -o gtdb-rs207.genomic-reps.dna.k31_not_in_sample.sig gtdb-rs207.genomic-reps.dna.k31.zip
# sketch the whole reference database in order to prep for a prefilter search
sourmash sketch dna -p k=31,scaled=1000,abund ../formatted_db.fasta -o formatted_db_merged.sig

# sourmash prefetch seems really slow for som reason, so let's use a gather and save the prefetch stuff, no need to worry about --threshold-bp since we won't actually use the gather results here, just the prefetch
# This will take a long time
sourmash gather --dna formatted_db_merged.sig gtdb-rs207.genomic-reps.dna.k31_not_in_sample.sig -o gather_formatted_db_merged_on_gtdb_not_in_sample.csv --save-prefetch-csv gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.csv --save-prefetch gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig
# the gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig file contains candidates from GTDB that are similar to the reference, but not in the sample. So now we need to:
# extract these candidate ones
sourmash sig extract --picklist gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.csv:match_md5:md5short -o gtdb-rs207.genomic-reps.dna.k31_not_in_sample_similar_to_merged_reference.sig gtdb-rs207.genomic-reps.dna.k31.zip
# then calculate the actual max ANI to the reference
python calculate_ANI_of_gtdb_to_reference.py

# then extract all the relevant sigs
cut -f2 in_gtdb_similar_to_EU_not_in_sample.tsv > in_gtdb_similar_to_EU_not_in_sample_md5.txt
sourmash sig split --output-dir sigs gather_formatted_db_merged_on_gtdb_not_in_sample_prefetch.sig --picklist in_gtdb_similar_to_EU_not_in_sample_md5.txt:gtdb_md5:md5


# get the sample sig
cp ../spike_in/36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig .

# get who the spikes are similar to
cp ../mut_thresh/in_gtdb_similar_to_EU_not_in_sample.tsv .
# sed squashes tabs
sed 's/\t/@/g' in_gtdb_similar_to_EU_not_in_sample.tsv > in_gtdb_similar_to_EU_not_in_sample_@.tsv

# Get the accessions of each of the spikes
sourmash sig collect -F csv -o sigs_manifest.csv sigs/*
cut -d',' -f2,3,10 sigs_manifest.csv | sed 's/"//g' | cut -d' ' -f1 | grep -v \# > sigs_md5_to_accession.txt

# find the location of the real GTDB genomes
cut -d',' -f3 sigs_md5_to_accession.txt | xargs -P 10 -I{} grep -m1 {} /data/shared_data/GTDB/gtdb_genomes_reps_r207/file_list.txt > sigs_gtdb_file_locations.txt

# add these to the accession list
paste -d',' <(cat sigs_md5_to_accession.txt) <(sed '1s/^/gtdb_file_location\n/' sigs_gtdb_file_locations.txt) > sigs_md5_to_accession_to_gtdb_location.txt

# For each of these GTDB genomes, get X coverage of the genome, sketch it, and then stick it in a folder
coverageValues=(".5" ".25" ".125" ".0625" ".03125" ".015625" ".0078125" ".00390625" ".001953125" ".0009765625")
for cov in ${coverageValues[@]}
do
	mkdir -p sigs_cov_${cov}
	mkdir -p sigs_cov_${cov}/reads
	mkdir -p spikes_cov_${cov}
	mkdir -p EU_on_spikes_cov_${cov}
done

# reduce the coverage
for cov in ${coverageValues[@]}
do
	for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location.txt`
	do
		md5short=$(echo ${line} | cut -d',' -f2)
		fileLoc=$(echo ${line} | cut -d',' -f4)
		# bbmap to coverage amount
		../../../KEGG_sketching_annotation/utils/bbmap/./randomreads.sh ref=${fileLoc} overwrite=t out=sigs_cov_${cov}/reads/${md5short}.fna coverage=0${cov}
	done
done
END
coverageValues=(".5" ".25" ".125" ".0625" ".03125" ".015625" ".0078125" ".00390625" ".001953125" ".0009765625")
# then sketch each one
for cov in ${coverageValues[@]}
do
        for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location.txt`
        do
                md5short=$(echo ${line} | cut -d',' -f2)
                fileLoc=$(echo ${line} | cut -d',' -f4)
		echo sourmash sketch dna -f -p k=31,scaled=1000,abund -o sigs_cov_${cov}/${md5short}.k=31.scaled=1000.DNA.dup=0.63.sig.zip sigs_cov_${cov}/reads/${md5short}.fna
	done | parallel -j 50
done

# then create all the spiked samples
for cov in ${coverageValues[@]}
do
        for line in `tail -n +2 sigs_md5_to_accession_to_gtdb_location.txt | head -n 10`
        do
                md5short=$(echo ${line} | cut -d',' -f2)
                fileLoc=$(echo ${line} | cut -d',' -f4)
                echo sourmash sig merge -o spikes_cov_${cov}/${md5short}_spiked_sample.sig.zip sigs_cov_${cov}/${md5short}.k=31.scaled=1000.DNA.dup=0.63.sig.zip 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig
        done | parallel -j 50
done

