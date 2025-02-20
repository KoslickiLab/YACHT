#!/bin/bash

#Create a folder for data files
mkdir data

### Download sample to data folder
fasterq-dump --concatenate-reads ERR5161194 -O data
fasterq-dump --concatenate-reads ERR10807863 -O data

### Download reference
datasets download genome accession GCF_016632365.1 --include genome --filename data/ncbi_dataset.zip
unzip data/ncbi_dataset.zip