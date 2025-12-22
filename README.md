# YACHT
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/KoslickiLab/YACHT/runTest.yml?logo=github&label=Python%20tests&link=https%3A%2F%2Fgithub.com%2FKoslickiLab%2FYACHT%2Factions)](https://github.com/KoslickiLab/YACHT/actions)
[![codecov](https://codecov.io/gh/KoslickiLab/YACHT/graph/badge.svg?token=AZD6LBFR5P)](https://codecov.io/gh/KoslickiLab/YACHT)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=KoslickiLab_YACHT&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=KoslickiLab_YACHT)
[![CodeQL](https://github.com/MichaelCurrin/badge-generator/workflows/CodeQL/badge.svg)](https://github.com/KoslickiLab/YACHT/actions?query=workflow%3ACodeQL "Code quality workflow status")
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/KoslickiLab/YACHT/blob/main/LICENSE.txt)

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on Average Nucleotide Identity (ANI). Identifying whether a specific microbe is actually present in a metagenomic sample is often complicated by sequencing noise, low-abundance organisms, and high genomic similarity between species. Traditional profiling tools rely on simple thresholds that can lead to high false-positive rates. Various cohorts can utilize YACHT: microbiome researchers dealing with low-biomass samples, synthetic biologists needing to validate the composition of mock communities, and genomics researchers identifying specific metagenome-assembled genomes (MAGs) of interest within vast sequencing datasets.

The associated publication can be found here: https://academic.oup.com/bioinformatics/article/40/2/btae047/7588873

And the preprint can be found at:  https://doi.org/10.1101/2023.04.18.537298.

Please cite via:

>Koslicki, D., White, S., Ma, C., & Novikov, A. (2024). YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. Bioinformatics, 40(2), btae047.

</br>

## Quick demonstration
We provide a demo to show how to use YACHT. Please follow the command lines below to try it out:

```bash
NUM_THREADS=64 # Adjust based on your machine's capabilities

cd demo # the 'demo' folder can be downloaded via command 'yacht download demo' if it doesn't exist

# build k-mer sketches for the query sample and ref genomes
yacht sketch sample --infile ./query_data/query_data.fq --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht sketch ref --infile ./ref_genomes --kmer 31 --scaled 1000 --outfile ref.sig.zip

# preprocess the reference genomes (training step)
yacht train --ref_file ref.sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./ --force

# run YACHT algorithm to check the presence of reference genomes in the query sample (inference step)
yacht run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads ${NUM_THREADS} --min_coverage_list 1 0.6 0.2 0.1 --outdir ./

# convert result to CAMI profile format (Optional)
yacht convert --yacht_output_dir ./results --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./
```

The output will be stored in the `results` folder containing:
- `result.xlsx`: An EXCEL file recording the presence of reference genomes with different spreadsheets given the minimum coverage of `1 0.6 0.2 0.1`.
- `result_all.txt`: A TXT file containing all unfiltered results for all user-given min_coverage values.

</br>

### Contents

- [YACHT](#yacht)
   * [Quick start](#quick-start)
   * [Installation](#installation)
      + [Conda Installation](#conda-installation)
      + [Manual installation](#manual-installation)
         - [Using Conda](#using-conda)
         - [Using Mamba](#using-mamba)
         - [Using Docker](#using-docker)
   * [Usage](#usage)
      + [YACHT Commands Overview](#yacht-commands-overview)
      + [YACHT workflow](#yacht-workflow)
      + [Creating sketches of your reference database genomes (yacht sketch ref)](#creating-sketches-of-your-reference-database-genomes-yacht-sketch-ref)
         - [Automatic download of reference sketches](#automatic-download-of-reference-sketches)
         - [Manual download of reference sketches](#manual-download-of-reference-sketches)
      + [Creating sketches of your sample (yacht sketch sample)](#creating-sketches-of-your-sample-yacht-sketch-sample)
      + [Preprocess the reference genomes (yacht train)](#preprocess-the-reference-genomes-yacht-train)
         - [Parameters](#parameters)
         - [Output](#output)
         - [Some pre-trained reference databases available on Zenodo  ](#some-pre-trained-reference-databases-available-on-zenodo)
      + [Run the YACHT algorithm (yacht run)](#run-the-yacht-algorithm-yacht-run)
         - [Parameters](#parameters-1)
         - [Output](#output-1)
      + [Convert YACHT result to other popular output formats (yacht convert)](#convert-yacht-result-to-other-popular-output-formats-yacht-convert)
         - [Parameters](#parameters-2)

## Installation

### Conda Installation

YACHT is available on Conda can be installed via the steps below to install：
```bash
# create conda environment
conda create -n yacht_env

# activiate environment
conda activate yacht_env

# install YACHT
conda install -c conda-forge -c bioconda yacht
```

### Manual installation
YACHT requires **Python >3.6** (and <3.12) with the following core genomics dependencies: `sourmash` (>=4.8.3), `sourmash_plugin_branchwater`, and `pytaxonkit`. The full list of dependencies can be found in the [environment configuration](https://github.com/KoslickiLab/YACHT/blob/main/env/yacht_env.yml). To ensure a clean and isolated workspace, we recommend using a virtual environment. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba), a faster alternative to Conda.
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
Using Dockerfile:
```
docker build --tag 'yacht' .
docker run -it --entrypoint=/bin/bash yacht -i
conda activate yacht_env

```
Using Act:

[Act](https://github.com/nektos/act). To run YACHT on docker, simply execute "act" from the main YACHT folder, or "act --container-architecture linux/amd64" if you are on MacOS system.

</br>

## Commands

YACHT can be run via the command line `yacht <module>`. The main modules include: `download`, `sketch`, `train`, `run`, and `convert`.

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


  + `pretrained_ref_db` can automatically download our pre-trained reference genome database that can be directly used as input for `yacht run` module.
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

- The `sketch` module (<ins>**note that it is a simple wrapper to `sourmash`**</ins>) has two submodules: `ref` and `sample`:
  
  + `ref` is used to sketch fasta files and make them as a reference database
  ```bash
  # Example for sketching multiple fasta files as reference genomes in a given folder
  yacht sketch ref --infile ./demo/ref_genomes --kmer 31 --scaled 1000 --outfile ref.sig.zip
  
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | infile            | the path to a input FASTQ file or a folder containing multiple FASTQ files |
  | kmer              | the length of k-mer |
  | scaled            | the scaled factor |
  | outfile           | the path to a output file |

  + `sample` is used to sketch the single-end or paired-end fasta file(s) and make it/them as a query sample.
  ```bash
  # Example for sketching a FASTA/Q file as a metagenomic example
  yacht sketch sample --infile ./query_data/query_data.fq --kmer 31 --scaled 1000 --outfile sample.sig.zip
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | infile            | the input FASTA/Q file(s). For paired-end reads, provide two files |
  | kmer              | the length of k-mer |
  | scaled            | the scaled factor |
  | outfile           | the path to a output file |

- The `train` module pre-reprocesses the given sketches of reference genomes (the `.zip` file) to identify and merge the "identical' genomes based on the given ANI threshold (e.g., --ani_threshold 0.95). For an example, please refer to the `yacht train` command in the "Quick start" section.

- The `run` module runs the YACHT algorithm to detect the presence of reference genomes in a given sample. For an example, please refer to the `yacht run` command in the "Quick start" section.

- The `convert` module can covert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn). For an example, please refer to the `yacht convert` command in the "Quick start" section.

## Workflow

This section introduces a brief workflow for using YACHT, summarized as:

1. **Create sketches of reference database genomes and samples:**
   
   `yacht sketch` samples compact representations of references or samples using `sourmash`.
2. **Preprocess the reference genomes:**
   
   `yacht train` preprocesses the reference genomes, merging those with high average nucleotide identity (ANI) into a single representative. 
3. **Run YACHT algorithm:** 
   
   `yacht run` executes the core YACHT algorithm to perform hypothesis testing and determine the presence or absence of organisms.
4. **Convert YACHT result to other output formats**
   
   `yacht convert` transforms the results into popular output formats like CAMI, BIOM, and GraphPhlAn.


### 1. Create sketches of reference database genomes and samples
#### Reference skeches


Use the command `yacht sketch` to generate sketches for both the samples and the reference genomes. Users must utilize sourmash to extract sketches from a reference database of microbial genomes. [sourmash Databases](https://sourmash.readthedocs.io/en/latest/databases.html) provide a variety of pre-formed databases of such sketches, or users can create a custom database using the `sourmash sketch` command on FASTA/FASTQ files of reference genomes (see the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/)). Other available databases include the [GTDB genomic representatives database](https://gtdb.ecogenomic.org/downloads). The sketches for samples must be generated using the same $k$-mer size and scale factor as those used for the reference database. The scale factor acts as an indicator of data compression, with smaller values being more appropriate for smaller datasets. 

We suggest trying with a pre-built reference sketches ([GTDB genomic representatives database](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip)):

```bash
yacht download default_ref_db --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder ./
```

Or
```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
```

For custom databases, you will need to create a Sourmash sketch Zipfile collection from the FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/)). Following commands accomplish it:

A single FASTA file with _one genome_ per record:

```bash
# This is equivalent to: sourmash sketch dna -f -p k=31,scaled=1000,abund --singleton <path to your multi-FASTA file> -o training_database.sig.zip
yacht sketch ref --infile <path to your multi-FASTA file> --kmer 31 --scaled 1000 --outfile training_database.sig.zip
```

A directory of FASTA files, one per genome:

```bash
# This is equivalent to: find <path of foler containg FASTA/FASTQ files> > dataset.csv; sourmash sketch fromfile dataset.csv -p dna,k=31,scaled=1000,abund -o training_database.sig.zip
yacht sketch ref --infile <path of foler containg FASTA/FASTQ files> --kmer 31 --scaled 1000 --outfile training_database.sig.zip
```

</br>


#### Sample skeches

This process should use the same k-mer size and scale factor that were used for the reference database. 

```bash
# For a single-end FASTA/Q file
# the command below is equivalent to: sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip <input FASTA/Q file>
yacht sketch sample --infile <input FASTA/Q file> --kmer 31 --scaled 1000 --outfile sample.sig.zip

# For pair-end FASTA/Q files, you need to separately specify two FASTA/Q files
# the command below is equivalent to: cat <FASTA/Q file 1> <FASTA/Q file 2> > combine.fastq (or combine.fasta); sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip combine.fastq (or combine.fasta)
yacht sketch sample --infile <FASTA/Q file 1> <FASTA/Q file 2> --kmer 31 --scaled 1000 --outfile sample.sig.zip
```

Note: Sourmash database offers three available k values (21, 31, and 51), allowing you to select the one that best suits your particular analytical needs. The scale factor serves as an indicator of data compression, and if your dataset is small, you might consider using a smaller value (corresponding to a higher portion of genomes retained in the sketch).


</br>

### 2. Preprocess the reference genomes (yacht train)

`yacht train` identifies and merges genomes that are roughly identical based on Average Nucleotide Identity (ANI). The module utilizes a fast algorithm written by C++ to preprocess the reference genomes. In our test with the GTDB representative genomes (r214) including `85,205` species-level genomes, YACHT takes around `12 minutes` and `52 GB` of RAM to preprocess them and generate the reference files on a Ubuntu 22.04.5 system using 64 threads. You can also use the pre-trained databases we built (see [here](#some-pre-trained-reference-databases-available-on-zenodo)) to skip this step.


```bash 
yacht train --ref_file gtdb-rs214-reps.k31.zip --ksize 31 --num_threads 64 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./
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

#### Some pre-trained reference databases available on Zenodo  

For convenience, we have provided some pre-trained reference database for the GenBank and GTDB genomes on [Zenodo](https://zenodo.org/communities/yacht?q=&l=list&p=1&s=10&sort=newest). If any of them is suitable for your study, you can simply run the following command to download it and skip the training step below. *Note*: download of pre-trained data is provided in the `yacht download` feature, please see [here](#yacht-commands-overview) for more details about `yacht download`.
```bash
# remember to replace <zendo_id> and <file_name> for your case before running it
curl --cookie zenodo-cookies.txt "https://zenodo.org/records/<zendo_id>/files/<file_name>?download=1" --output <file_name>

# Example
# curl --cookie zenodo-cookies.txt "https://zenodo.org/records/10113534/files/genbank-2022.03-archaea-k31_0.80_pretrained.zip?download=1" --output genbank-2022.03-archaea-k31_0.80_pretrained.zip
```


</br>

### 3. Run the YACHT algorithm (yacht run)

After this, you are ready to perform the hypothesis test via `yacht run` for each organism in your reference database. This can be accomplished with something like:

```bash
yacht run --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --num_threads 64 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --outdir ./
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
| --outdir                          | path to output location where the `results` folder will be created (default: current working directory) |

#### Output

The output will be stored in the `results` folder at the specified `--outdir` location, containing:

| File                  | Content                                                      |
| --------------------- | ------------------------------------------------------------ |
| result.xlsx           | An EXCEL file with filtered results for each min_coverage value (one sheet per value) |
| result_all.txt   | A TXT file containing all unfiltered results for all user-given min_coverage values |

The column descriptions can be found [here](docs/column_descriptions.csv). The most important are the following:

* `organism_name`: The name of the organism
* `in_sample_est`: A boolean value either False or True: if False, there was not enough evidence to claim this organism is present in the sample. 
* `p_vals`: Probability of observing this or more extreme result at the given ANI threshold, assuming the null hypothesis.

Other interesting columns include:

* `num_exclusive_kmers_to_genome`: How many k-mers were found in this organism and no others
* `num_matches`: How many k-mers were found in this organism and the sample
* `acceptance_threshold_*`: How many k-mers must be found in this organism to be considered "present" at the given ANI threshold. Hence, `in_sample_est` is True if `num_matches` >= `acceptance_threshold_*` (adjusting by coverage if desired).
* `alt_confidence_mut_rate_*`: What the mutation rate (1-ANI) would need to be to get your false positive to match the false negative rate of 1-`significance` (adjusting by coverage if desired).

</br>

### 4. Convert YACHT result to other popular output formats (yacht convert)

When we get the results folder from `yacht run`, you can run `yacht convert` to covert the YACHT result to other popular output formats (Currently, only `cami`, `biom`, `graphplan` are supported).

__Note__: Before you run `yacht convert`, you need to prepare a TSV file `genome_to_taxid.tsv` containing two columns: genome ID (genome_id) and its corresponding taxid (taxid). An example can be found [here](demo/toy_genome_to_taxid.tsv). You need to prepare it according to the reference database genomes you used. 

Then you are ready to run `yacht convert` with something like:
```bash
yacht convert --yacht_output_dir './results' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir ./
```

#### Parameters

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --yacht_output_dir    | the path to the `results` folder generated by `yacht run` (containing `result.xlsx`) |
| --sheet_name      | specify which spreadsheet result you want to covert from     |
| --genome_to_taxid | the path to the location of `genome_to_taxid.tsv` you prepared |
| --mode            | specify to which output format you want to convert (e.g., 'cami', 'biom', 'graphplan')
| --sample_name     | A random name you would like to show in header of the cami file. Default: Sample1.' |
| --outfile_prefix  | the prefix of the output file. Default: result | 
| --outdir          | the path to output directory where the results will be genreated |



