#!/usr/bin/env bash
set -e
set -u
set -o pipefail

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

# then create all the spiked samples
for file in `ls sigs/*`
do
	name=$(basename ${file})
	echo sourmash sig merge -o spikes/${name}_spike.sig ${file} 36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig
done | parallel -j 100


