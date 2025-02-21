#!/bin/bash

#Create a folder for data files
mkdir data

### Download sample to data folder
fasterq-dump --concatenate-reads SRR32008482  -O data
fasterq-dump --concatenate-reads SRR32008483  -O data

### Download MAG of interest
datasets download genome accession GCA_017506175.1 --include genome --filename data/.
cd data
unzip ncbi_dataset
cd ../
