#!/bin/bash
set -e
set -u
set -o pipefail
i=0
while read name
do
	i=$((i+1))
	echo sourmash sig extract -f --include-db-pattern "\"${name}\"" /data/shared_data/sourmash_data/wxf9z/googledrive/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip -o sigs/${i}.sig
done < in_gtdb_similar_to_EU_not_in_sample_names.txt
