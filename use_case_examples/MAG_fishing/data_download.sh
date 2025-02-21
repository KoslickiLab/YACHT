#!/bin/bash

#Create a folder for data files
mkdir data

### Download sample to data folder
fasterq-dump --concatenate-reads SRR32008482  -O data
fasterq-dump --concatenate-reads SRR32008483  -O data

### Download pre-reference signatures for k=31 and k=51
#k-size=31
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.zip --directory-prefix=data
#k-size=51
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k51.zip --directory-prefix=data

### Download MAG of interest
datasets download genome accession GCA_017506175.1 --include genome --filename data/.
cd data
unzip ncbi_dataset
cd ../
