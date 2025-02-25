#!/bin/bash

#Create a folder for data files
mkdir data

### Download sample to data folder
fasterq-dump --concatenate-reads SRR32008482  -O data

### Download pre-reference signatures for k=21,k=31,and k=51
#k-size=21
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k21.zip --directory-prefix=data
#k-size=31
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.zip --directory-prefix=data
#k-size=51
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k51.zip --directory-prefix=data
