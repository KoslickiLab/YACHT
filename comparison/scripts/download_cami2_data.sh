#! /bin/bash
# Download CAMI2 data from https://data.cami-challenge.org/participate
# Usage: download_cami2_data.sh <output_dir> <cpu_num>

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: download_cami2_data.sh <output_dir> <cpu_num>"
    exit 1
fi
output_dir=$1
cpu_num=$2

# create a CAMI_data directory for the output
if [ ! -d $output_dir/CAMI_data ]; then
    mkdir -p $output_dir/CAMI_data
fi

# download the data seperately for Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset
# Rhizosphere challenge
if [ ! -d $output_dir/CAMI_data/rhizosphere_data ]; then
    mkdir -p $output_dir/CAMI_data/Rhizosphere_data;
    parallel -j $cpu_num wget -P $output_dir/CAMI_data/rhizosphere_data https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_sample_{}_reads.tar.gz ::: {0..20};
    ## download the OPAL result
    wget -P $output_dir/CAMI_data/rhizosphere_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/rhizosphere_dataset/results/OPAL_short_long/results.tsv;
    ## download the ground truth profile
    wget -P $output_dir/CAMI_data/rhizosphere_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/rhizosphere_dataset/data/ground_truth/gs_rhizosphere.profile;
fi


# Clinical pathogen detection challenge
if [ ! -d $output_dir/CAMI_data/pathogen_detection_data ]; then
    mkdir -p $output_dir/CAMI_data/pathogen_detection_data;
    wget -P $output_dir/CAMI_data/pathogen_detection_data https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz;
    tar zxvf $output_dir/CAMI_data/pathogen_detection_data/patmgCAMI2.tar.gz -C $output_dir/CAMI_data/pathogen_detection_data;
    mv $output_dir/CAMI_data/pathogen_detection_data/patmg_CAMI2/* $output_dir/CAMI_data/pathogen_detection_data;
    rm -r $output_dir/CAMI_data/pathogen_detection_data/patmg_CAMI2;
    rm $output_dir/CAMI_data/pathogen_detection_data/patmgCAMI2.tar.gz;
fi


# Challenge Marine Dataset
if [ ! -d $output_dir/CAMI_data/marine_data ]; then
    mkdir -p $output_dir/CAMI_data/marine_data
    parallel -j $cpu_num wget -P $output_dir/CAMI_data/marine_data https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_{}_reads.tar.gz ::: {0..9};
    # download the OPAL result
    wget -P $output_dir/CAMI_data/marine_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/results/OPAL_default_short_read_samples/results.tsv;
    # download the ground truth profile
    wget -P $output_dir/CAMI_data/marine_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/data/ground_truth/gs_marine_short.profile;
fi

# Strain Madness Dataset
if [ ! -d $output_dir/CAMI_data/strain_madness_data ]; then
    mkdir -p $output_dir/CAMI_data/strain_madness_data;
    parallel -j $cpu_num wget -P $output_dir/CAMI_data/strain_madness_data https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_sample_{}_reads.tar.gz ::: {0..99};
    # download the OPAL result
    wget -P $output_dir/CAMI_data/strain_madness_data/opal_results https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/strain_madness_dataset/results/OPAL_default_short_read_samples/results.tsv;
    # download the ground truth profile
    wget -P $output_dir/CAMI_data/strain_madness_data/ground_truth https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/strain_madness_dataset/data/ground_truth/gs_strain_madness_short_long.profile;
fi