#!/bin/bash
set -e
set -u
set -o pipefail

#This script will take all of the GTDB organisms that are similar to the EU training matrix, but not in the sample, then spike them into the sample, saving the *.sig file for later `sourmash gather` usage and parsing

spikeDir="spikes"
sampleSketch="36116.SZAXPI030664-33.clean.trim.rmhost.1.fq.sig"
gtdbRef="prefetch_formatted_db_to_gtdb4_not_in_sample.sig"
# for each entry, get the signature and merge it with the sample
for file in `ls sigs/*`
do
	name=$(basename ${file})
	echo sourmash sig merge -o ${spikeDir}/${name}_spike.sig ${file} ${sampleSketch}
done | parallel -j 200
