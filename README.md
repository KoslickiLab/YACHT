# YACHT
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/KoslickiLab/YACHT/runTest.yml?logo=github&label=Python%20tests&link=https%3A%2F%2Fgithub.com%2FKoslickiLab%2FYACHT%2Factions)](https://github.com/KoslickiLab/YACHT/actions)
[![codecov](https://codecov.io/gh/KoslickiLab/YACHT/graph/badge.svg?token=AZD6LBFR5P)](https://codecov.io/gh/KoslickiLab/YACHT)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=KoslickiLab_YACHT&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=KoslickiLab_YACHT)
[![CodeQL](https://github.com/MichaelCurrin/badge-generator/workflows/CodeQL/badge.svg)](https://github.com/KoslickiLab/YACHT/actions?query=workflow%3ACodeQL "Code quality workflow status")
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/KoslickiLab/YACHT/blob/main/LICENSE.txt)

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).

The associated preprint can be found at:  https://doi.org/10.1101/2023.04.18.537298. Please cite via:

>Koslicki, D., White, S., Ma, C., & Novikov, A. (2023). YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. bioRxiv, 2023-04.

</br>

## Quick start
We provide a demo to show how to use YACHT. Please follow the command lines below to try it out:

```bash
NUM_THREADS=64 # Adjust based on your machine's capabilities

cd demo # the 'demo' folder can be downloaded via command 'yacht download demo' if it doesn't exist

# build k-mer sketches for the query sample and ref genomes
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip query_data/query_data.fq
sourmash sketch fromfile ref_paths.csv -p dna,k=31,scaled=1000,abund -o ref.sig.zip --force-output-already-exists

# preprocess the reference genomes (training step)
yacht train --ref_file ref.sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./ --force

# run YACHT algorithm to check the presence of reference genomes in the query sample (inference step)
yacht run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads ${NUM_THREADS} --min_coverage_list 1 0.6 0.2 0.1 --out ./result.xlsx

# convert result to CAMI profile format (Optional)
yacht convert --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./
```

There will be an output EXCEL file `result.xlsx` recoding the presence of reference genomes with different spreadsheets given the minimum coverage of `1 0.6 0.2 0.1`.

</br>

### Contents

- [Installation](#installation)
  * [Conda Installation](#conda-installation)
  * [Manual installation](#manual-installation)
    + [Using Conda](#using-conda)
    + [Using Mamba](#using-mamba)
    + [Using Docker](#using-docker)
- [Usage](#usage)
  * [YACHT Commands Overview](#yacht-commands-overview)
  * [YACHT workflow](#yacht-workflow)
  * [Creating sketches of your reference database genomes](#creating-sketches-of-your-reference-database-genomes)
    + [Automatic download of reference sketches](#automatic-download-of-reference-sketches)
    + [Manual download of reference sketches](#manual-download-of-reference-sketches)
  * [Creating sketches of your sample](#creating-sketches-of-your-sample)
    + [Parameters](#parameters)
    + [Output](#output)
  * [Preprocess the reference genomes](#preprocess-the-reference-genomes-yacht-train)
    + [Parameters](#parameters-1)
    + [Output](#output-1)
    + [Some pre-trained reference databases available on Zenodo](#some-pre-trained-reference-databases-available-on-zenodo)
  * [Run the YACHT algorithm](#run-the-yacht-algorithm)
    + [Parameters](#parameters-2)
    + [Output](#output-2)
  * [Convert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn)](#convert-yacht-result-to-other-popular-output-formats-eg-cami-profiling-format-biom-format-graphplan)
    + [Parameters](#parameters-3)

## Installation

**Please note YACHT does not currently support MacOS. However, we are actively working on developing compatibility for this operating system and hope to have it available soon. During this time, we provide a docker container (see `using docker` section below) for those who need to run YACHT on MacOS.**

### Conda Installation

A Conda package for YACHT will be available soon. Once it is available, YACHT can be installed via the steps below：
```bash
# create conda environment
conda create -n yacht_env

# activiate environment
conda activate yacht_env

# install YACHT
conda install -c bioconda yacht
```

### Manual installation
YACHT requires Python 3.6 or higher. We recommend using a virtual environment to ensure a clean and isolated workspace. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

#### Using Conda
To create your Conda environment and install YACHT, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/YACHT.git
cd YACHT

# Create a new virtual environment named 'yacht_env'
conda env create -f env/yacht_env.yml

# Activate the newly created environment
conda activate yacht_env

# Install YACHT within the environment
pip install .
```

#### Using Mamba
If you prefer using Mamba instead of Conda, just simply repalce `conda` with `mamba` in the above commands.

#### Using Docker
If you prefer running YACHT on MacOS, you can choose to use docker with [Act](https://github.com/nektos/act). To run YACHT on docker, simply execute "act" from the main YACHT folder, or "act --container-architecture linux/amd64" if you are on MacOS system.

</br>

## Usage

### YACHT Commands Overview
YACHT can be run via the command line `yacht <module>`. Now it has three four main modules: `download`, `train`, `run`, and `convert`.

- The `download` module has three submodules: `demo`, `default_ref_db`, and `pretrained_ref_db`:
  
  + `demo` can automatically download the demo files to a specified folder:
  ```bash
  # Example
  yacht download demo --outfolder ./demo
  ```
  + `default_ref_db` can automatically download pre-generated sketches of reference genomes from GTDB or GenBank as our input reference databases.
  ```bash
  # Example for downloading the k31 sketches of representative genomes of GTDB rs214 version 
  yacht download default_ref_db --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder ./
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | database          | two options for default reference databases: 'genbank' or 'gtdb' |
  | db_version        | the version of database, options: "genbank-2022.03", "rs202", "rs207", "rs214" |
  | ncbi_organism     | the NCBI organism for the NCBI reference genome, options: "archaea", "bacteria", "fungi", "virus", "protozoa"|
  | gtdb_type         | for GTDB database, chooses "representative" genome version or "full" genome version |
  | k                 | the length of k-mer |
  | outfolder         | the path to a folder where the downloaded file is expected to locate |


  + `pretrained_ref_db` can automatically download our pre-trained reference genome database that can be directly used as input for `yacht train` module.
  ```bash
  # Example for downloading the pretrained reference database that was trained from GTDB rs214 representative genomes with k=31 and ani_threshold=0.9995
  yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 31 --ani_thresh 0.9995 --outfolder ./
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | database          | two options for default reference databases: 'genbank' or 'gtdb' |
  | db_version        | the version of database, options: "genbank-2022.03", "rs214" |
  | ncbi_organism     | the NCBI organism for the NCBI reference genome, options: "archaea", "bacteria", "fungi", "virus", "protozoa"|
  | ani_thresh      | the cutoff by which two organisms are considered indistinguishable (default: 0.95) |
  | k                 | the length of k-mer |
  
  | outfolder         | the path to a folder where the downloaded file is expected to locate |

- The `train` module pre-reprocesses the given sketches of reference genomes (the `.zip` file) to identify and merge the "identical' genomes based on the given ANI threshold (e.g., --ani_threshold 0.95). For an example, please refer to the `yacht train` command in the "Quick start" section.

- The `run` module runs the YACHT algorithm to detect the presence of reference genomes in a given sample. For an example, please refer to the `yacht run` command in the "Quick start" section.

- The `convert` module can covert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn). For an example, please refer to the `yacht convert` command in the "Quick start" section.

### YACHT workflow

This section simply introduces the analysis workflow for YACHT:

1. **Create Sketches of Your Reference Database Genomes and Your Sample:**
   - This involves using [sourmash](https://sourmash.readthedocs.io/en/latest/) to generate compact representations (sketches) of genomic data for efficient comparison and analysis.
2. **Preprocess the Reference Genomes:**
   - This is the training step of YACHT, aiming to identify and merge the "identical" genomes based on Average Nucleotide Identity (`ANI`) using the `ani_thresh` parameter. 
3. **Run YACHT algorithm:** 
   - This step involves running the YACHT algorithm to detect the presence of reference genomes in your sample.
4. **Convert YACHT result to other output formats**
   - This step is optional if you prefer other output formats (e.g., CAMI profiling format, BIOM format) for the downstream analysis.


For each step of this workflow, please see more detailed description in the sections below.

</br>

### Creating sketches of your reference database genomes

You will need a reference database in the form of sourmash sketches of a collection of microbial genomes. There are a variety of pre-created databases available at: https://sourmash.readthedocs.io/en/latest/databases.html. Our code uses the "Zipfile collection" format, and we suggest using the [GTDB genomic representatives database](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip):

#### Automatic download of reference sketches
```bash
yacht download default_ref_db --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder ./
```

#### Manual download of reference sketches
```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
```

If you want to use a custom database, you will need to create a Sourmash sketch Zipfile collection from the FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/) for details). In brief, this can be accomplished via the following commands:

If you have a single FASTA file with _one genome_ per record:

```bash
sourmash sketch dna -f -p k=31,scaled=1000,abund --singleton <your multi-FASTA file> -o training_database.sig.zip
```

If you have a directory of FASTA files, one per genome:

```bash
## Method 1 (suggested)
# put all full paths of FASTA/FASTQ file into a file, one path per line
find <path of foler containg FASTA/FASTQ files> > dataset.csv
sourmash sketch fromfile dataset.csv -p dna,k=31,scaled=1000,abund -o ../training_database.sig.zip
# cd back to YACHT

## Method 2
# cd into the relevant directory
sourmash sketch dna -f -p k=31,scaled=1000,abund *.fasta -o ../training_database.sig.zip
# cd back to YACHT
```

</br>


### Creating sketches of your sample

Creating a sketch of your sample metagenome is an essential step in the YACHT workflow. This process involves using the same k-mer size and scale factor that were used for the reference database. You can use the following commands to implement this step:

```bash
# For a single-end FASTA/Q file
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip <input FASTA/Q file>

# For pair-end FASTA/Q files, you need to combine them into a single file first
cat <FASTA/Q file 1> <FASTA/Q file 2> > combine.fastq (or combine.fasta)
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip combine.fastq (or combine.fasta)
```

#### Parameters

Sourmash database offers three available k values (21, 31, and 51), allowing you to select the one that best suits your particular analytical needs. The scale factor serves as an indicator of data compression, and if your dataset is small, you might consider using a smaller value (corresponding to a higher portion of genomes retained in the sketch).

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| k=31              | the length of k-mer                                          |
| scaled=1000       | fraction of k-mers to be kept in the sketch                  |
| -o sample.sig.zip | arbitrary output name for sketch files (must be in zip format) |

#### Output

In the two preceding steps, you will obtain a k-mer sketch file in zip format (in each step) for the reference data and the input sample.

| File                      | Content                                                     |
| ------------------------- | ----------------------------------------------------------- |
| gtdb-rs214-reps.k31.zip   | Pre-built k-mer sketch file for GTDB representative genomes |
| training_database.sig.zip | K-mer sketch file for your selected reference genomes       |
| sample.sig.zip            | K-mer sketch file for your input sample                     |


</br>

### Preprocess the reference genomes (yacht train)

**Warning: the training process is time-consuming on large database**

In our benchmark with `GTDB representive genomes`, it takes `15 minutes` using `16 threads, 50GB of MEM` on a system equipped with a `3.5GHz AMD EPYC 7763 64-Core Processor`. You can use the pre-trained database (see [here](#some-pre-trained-reference-databases-available-on-zenodo)) to skip this step. The processing time can be significant when executed on GTDB all genomes OR with limited resources. If only part of genomes are needed, one may use `sourmash sig` command to extract signatures of interests only. 

</br>

The command `yacht train` extracts the sketches from the Zipfile-format reference database, and then turns them into a form usable by YACHT. In particular, it removes one of any two organisms that have ANI greater than the user-specified threshold as these two organisms are too close to be "distinguishable".

```bash 
yacht train --ref_file gtdb-rs214-reps.k31.zip --ksize 31 --num_threads 32 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./
```

#### Parameters

The most important parameter of this script is `--ani_thresh`: this is average nucleotide identity (ANI) value equal to or below which two organisms are considered distinct. For example, if `--ani_thresh` is set to 0.95, then two organisms with ANI > 0.95 will be considered indistinguishable. For the organisms with ANI > 0.95, only the one with the largest number of unique kmers will be kept. If there is a tie in the number of unique kmers, one organism will be randomly selected. The default value of `--ani_thresh` is 0.95. The `--ani_thresh` value chosen here must match the one chosen for the YACHT algorithm (see below).  

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --ref_file        | the path to the sourmash signature database zip file |
| --ksize           | the length of k-mer, must match the k size used in previous sketching steps (default: 31) |
| --num_threads     | the number of threads to use for parallelization (default: 16) |
| --ani_thresh      | the cutoff by which two organisms are considered indistinguishable (default: 0.95) |
| --prefix          | the prefix for output folders and files (see details below) |
| --outdir          | the path to output directory where the results and intermediate files will be genreated |

#### Output

| File (names starting with prefix)     | Content                                                      |
| ------------------------------------- | ------------------------------------------------------------ |
| _config.json                          | A JSON file stores the required information needed to run the next YACHT algorithm |
| _manifest.tsv                         | A TSV file contains organisms and their relevant info after removing the similar ones |
| _removed_orgs_to_corr_orgas_mapping.tsv   | A TSV file with two columns: removed organism names ('removed_org') and their similar genomes ('corr_orgs')| 

#### Some pre-trained reference databases available on Zenodo  

For convenience, we have provided some pre-trained reference database for the GenBank and GTDB genomes on [Zenodo](https://zenodo.org/communities/yacht?q=&l=list&p=1&s=10&sort=newest). If any of them is suitable for your study, you can simply run the following command to download it and skip the training step below:
```bash
# remember to replace <zendo_id> and <file_name> for your case before running it
curl --cookie zenodo-cookies.txt "https://zenodo.org/records/<zendo_id>/files/<file_name>?download=1" --output <file_name>

# Example
# curl --cookie zenodo-cookies.txt "https://zenodo.org/records/10113534/files/genbank-2022.03-archaea-k31_0.80_pretrained.zip?download=1" --output genbank-2022.03-archaea-k31_0.80_pretrained.zip
```

**Please note that if you plan to use these pre-trained reference databases, once you download and unzip it. You need to change the paths within the config json file (e.g., gtdb-rs214-reps.k31_0.9995_config.json) to the correct paths in your machine.**

</br>

### Run the YACHT algorithm (yacht run)

After this, you are ready to perform the hypothesis test via `yacht run` for each organism in your reference database. This can be accomplished with something like:

```bash
yacht run --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --num_threads 32 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result.xlsx
```

#### Parameters

The `--significance` parameter is basically akin to your confidence level: how sure do you want to be that the organism is present? Higher leads to more false negatives, lower leads to more false positives. 

The `--min_coverage_list` parameter dictates a list of `min_coverage` which indicates what percentage (value in `[0,1]`) of the distinct k-mers (think: whole genome) must have been sequenced and present in my sample to qualify as that organism as being "present." Setting this to 1 is usually safe, but if you have a very low coverage sample, you may want to lower this value. Setting it higher will lead to more false negatives, setting it lower will lead to more false positives (pretty rapidly).

| Parameter                               | Explanation                                                  |
| --------------------------------------- | ------------------------------------------------------------ |
| --json      | the path to a json file generated by the `make_training_data_from_sketches.py` script (see above) |
| --significance                          | minimum probability of individual true negative (default: 0.99) |
| --num_threads                           | the number of threads to use for parallelization (default: 16) |
| --keep_raw                              | keep the raw result (i.e. `min_coverage=1`) no matter if the user specifies it |
| --show_all                              | Show all organisms (no matter if present) |
| --min_coverage_list                     | a list of `min_coverage` values, see more detailed description above (default: 1, 0.5, 0.1, 0.05, 0.01) |
| --out                          | path to output excel result (default: './result.xlsx') |

#### Output

The output file will be an EXCEL file; column descriptions can be found [here](docs/column_descriptions.csv). The most important are the following:

* `organism_name`: The name of the organism
* `in_sample_est`: A boolean value either False or True: if False, there was not enough evidence to claim this organism is present in the sample. 
* `p_vals`: Probability of observing this or more extreme result at the given ANI threshold, assuming the null hypothesis.

Other interesting columns include:

* `num_exclusive_kmers_to_genome`: How many k-mers were found in this organism and no others
* `num_matches`: How many k-mers were found in this organism and the sample
* `acceptance_threshold_*`: How many k-mers must be found in this organism to be considered "present" at the given ANI threshold. Hence, `in_sample_est` is True if `num_matches` >= `acceptance_threshold_*` (adjusting by coverage if desired).
* `alt_confidence_mut_rate_*`: What the mutation rate (1-ANI) would need to be to get your false positive to match the false negative rate of 1-`significance` (adjusting by coverage if desired).

</br>

### Convert YACHT result to other popular output formats (yacht convert)

When we get the EXCEL result file from run_YACHT.py, you can run `yacht convert` to covert the YACHT result to other popular output formats (Currently, only `cami`, `biom`, `graphplan` are supported).

__Note__: Before you run `yacht convert`, you need to prepare a TSV file `genome_to_taxid.tsv` containing two columns: genome ID (genome_id) and its corresponding taxid (taxid). An example can be found [here](demo/toy_genome_to_taxid.tsv). You need to prepare it according to the reference database genomes you used. 

Then you are ready to run `yacht convert` with something like:
```bash
yacht convert --yacht_output 'result.xlsx' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir ./
```

#### Parameters

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --yacht_output    | the path to the output excel file generated by `run_YACHT.py` |
| --sheet_name      | specify which spreadsheet result you want to covert from     |
| --genome_to_taxid | the path to the location of `genome_to_taxid.tsv` you prepared |
| --mode            | specify to which output format you want to convert (e.g., 'cami', 'biom', 'graphplan')
| --sample_name     | A random name you would like to show in header of the cami file. Default: Sample1.' |
| --outfile_prefix  | the prefix of the output file. Default: result | 
| --outdir          | the path to output directory where the results will be genreated |


