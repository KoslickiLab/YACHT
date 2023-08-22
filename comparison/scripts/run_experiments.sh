#! /bin/bash
# run benchmarking experiments based on CAMI2 data for YACHT algorithm
# Usage:    

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run_experiments.sh <yacht_repo_loc> <output_dir> <cpu_num>"
    exit 1
fi
yacht_repo_loc=$1
output_dir=$2
cpu_num=$3

## Run YACHT on CAMI2 data (e.g., Rhizosphere challenge, Clinical pathogen detection challenge, Challenge Marine Dataset, Strain Madness Dataset)
# Download GTDB genomics representative database for YACHT
cd $yacht_repo_loc
if [ ! -d $yacht_repo_loc/gtdb_reference ]; then
    mkdir -p $yacht_repo_loc/gtdb_reference;
    wget -P $yacht_repo_loc/gtdb_reference https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
fi

# Crate a reference dictionary matrix for YACHT
python $yacht_repo_loc/make_training_data_from_sketches.py --ref_file $yacht_repo_loc/gtdb_reference/gtdb-rs214-reps.k31.zip --ksize 31 --out_prefix 'gtdb_ani_thresh_0.95' --ani_thresh 0.95

## run YACHT on CAMI2 data
# create a "results" directory for the output
if [ ! -d $output_dir/results ]; then
    mkdir -p $output_dir/results
fi
if [ ! -d $output_dir/results/YACHT_results ]; then
    mkdir -p $output_dir/results/YACHT_results
fi

# run YACHAT on Rhizosphere challenge data
if [ ! -d $output_dir/results/YACHT_results/rhizosphere_data ]; then
    mkdir -p $output_dir/results/YACHT_results/rhizosphere_data;
    ## create sketches of samples
    if [ ! -d $output_dir/CAMI_data/rhizosphere_data/sketches ]; then
        mkdir -p $output_dir/CAMI_data/rhizosphere_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund --outdir $output_dir/CAMI_data/rhizosphere_data/sketches -o rhimgCAMI2_sample_{}.sig.zip $output_dir/CAMI_data/rhizosphere_data/rhimgCAMI2_sample_{}_reads.tar.gz ::: {0..20};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $output_dir/CAMI_data/rhizosphere_data/sketches/rhimgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $output_dir/results/YACHT_results/rhizosphere_data/rhimgCAMI2_sample_{}.csv ::: {0..20};
fi

# run YACHT on Clinical pathogen detection challenge data
if [ ! -d $output_dir/results/YACHT_results/pathogen_detection_data ]; then
    mkdir -p $output_dir/results/YACHT_results/pathogen_detection_data;
    ## create sketches of samples
    if [ ! -d $output_dir/CAMI_data/pathogen_detection_data/sketches ]; then
        mkdir -p $output_dir/CAMI_data/pathogen_detection_data/sketches;
        sourmash sketch dna -f -p k=31,scaled=1000,abund --outdir $output_dir/CAMI_data/pathogen_detection_data/sketches -o patmgCAMI2.sig.zip $output_dir/CAMI_data/pathogen_detection_data/*.tar.gz;
    fi
    ## run YACHT on the samples
    python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $output_dir/CAMI_data/pathogen_detection_data/sketches/patmgCAMI2.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $output_dir/results/YACHT_results/pathogen_detection_data/patmgCAMI2.csv;
fi


# run YACHT on Challenge Marine Dataset
if [ ! -d $output_dir/results/YACHT_results/marine_data ]; then
    mkdir -p $output_dir/results/YACHT_results/marine_data;
    ## create sketches of samples
    if [ ! -d $output_dir/CAMI_data/marine_data/sketches ]; then
        mkdir -p $output_dir/CAMI_data/marine_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund --outdir $output_dir/CAMI_data/marine_data/sketches -o marmgCAMI2_sample_{}.sig.zip $output_dir/CAMI_data/marine_data/marmgCAMI2_sample_{}_reads.tar.gz ::: {0..9};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $output_dir/CAMI_data/marine_data/sketches/marmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $output_dir/results/YACHT_results/marine_data/marmgCAMI2_sample_{}.csv ::: {0..9};
fi


# run YACHT on Strain Madness Dataset
if [ ! -d $output_dir/results/YACHT_results/strain_madness_data ]; then
    mkdir -p $output_dir/results/YACHT_results/strain_madness_data;
    ## create sketches of samples
    if [ ! -d $output_dir/CAMI_data/strain_madness_data/sketches ]; then
        mkdir -p $output_dir/CAMI_data/strain_madness_data/sketches;
        parallel -j $cpu_num sourmash sketch dna -f -p k=31,scaled=1000,abund --outdir $output_dir/CAMI_data/strain_madness_data/sketches -o strmgCAMI2_sample_{}.sig.zip $output_dir/CAMI_data/strain_madness_data/strmgCAMI2_sample_{}_reads.tar.gz ::: {0..99};
    fi
    ## run YACHT on the samples
    parallel -j $cpu_num python $yacht_repo_loc/run_YACHT.py --ref_matrix $yacht_repo_loc/gtdb_ani_thresh_0.95_ref_matrix_processed.npz --sample_file $output_dir/CAMI_data/strain_madness_data/sketches/strmgCAMI2_sample_{}.sig.zip --ksize 31 --ani_thresh 0.95 --significance 0.99 --min_coverage 1 --outfile $output_dir/results/YACHT_results/strain_madness_data/strmgCAMI2_sample_{}.csv ::: {0..99};
fi

