# YACHT
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/KoslickiLab/YACHT/runTest.yml?logo=github&label=Python%20tests&link=https%3A%2F%2Fgithub.com%2FKoslickiLab%2FYACHT%2Factions)](https://github.com/KoslickiLab/YACHT/actions)
[![codecov](https://codecov.io/gh/KoslickiLab/YACHT/graph/badge.svg?token=AZD6LBFR5P)](https://codecov.io/gh/KoslickiLab/YACHT)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=KoslickiLab_YACHT&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=KoslickiLab_YACHT)
[![CodeQL](https://github.com/MichaelCurrin/badge-generator/workflows/CodeQL/badge.svg)](https://github.com/KoslickiLab/YACHT/actions?query=workflow%3ACodeQL "Code quality workflow status")
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/KoslickiLab/YACHT/blob/main/LICENSE.txt)

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI). Now featuring enhanced capabilities for downloading demo data, genome sketches, and pretrained models, and for streamlining the sketch creation and result standardization processes.

The associated preprint can be found at:  https://doi.org/10.1101/2023.04.18.537298. Please cite via:

>Koslicki, D., White, S., Ma, C., & Novikov, A. (2023). YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. bioRxiv, 2023-04.



## Quick start
Get started with YACHT using our demo:

```bash
NUM_THREADS=16  # Adjust based on your machine's capabilities

# Download demo data
yacht_download_demo --output demo

cd demo

# Create sketches for the query sample and reference genomes
yacht_sketch_sample --infile query_data/query_data.fq --k 31 --scaled 1000 --outfile sample.sig.zip
yacht_sketch_genomes --infile ref_genomes --k 31 --scaled 1000 --outfile ref.sig.zip

# Preprocess reference genomes (training step)
yacht_training --ref_file ref.sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./

# Run YACHT
yacht_run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads ${NUM_THREADS} --min_coverage_list 1 0.6 0.2 0.1 --out result.xlsx

# Convert result to other formats (Optional)
yacht_standardize --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./

```

There will be an output EXCEL file `result.xlsx` recoding the presence of reference genomes with different spreadsheets given the minimum coverage of `1 0.6 0.2 0.1`.


## Installation

### Conda Installation

A Conda package for YACHT will be available soon. Until then, YACHT can be installed manually following the instructions below.

### Manual Installation

YACHT requires Python 3 or higher. We recommend using a virtual environment to ensure a clean and isolated workspace. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/conda-forge/miniforge), a faster alternative to Conda.

#### Using Mamba

To create and manage your virtual environment using Mamba, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/YACHT.git
cd YACHT

# Create a new virtual environment named 'yacht_env'
mamba env create -f env/yacht_env.yml

# Activate the newly created environment
mamba activate yacht_env

# Install YACHT within the environment
pip install .
```

#### Using Conda

If you prefer using Conda over Mamba, simply replace `mamba` with `conda` in the above commands.

## Usage

The YACHT workflow encompasses the following steps:

1. **Creating Sketches of Your Reference Database Genomes and Your Sample:**
   - This involves generating compact representations (sketches) of genomic data for efficient comparison and analysis.

2. **Preprocessing the Reference Genomes:**
   - This step involves removing genomes that are too similar based on Average Nucleotide Identity (ANI), using the `ani_thresh` parameter. This ensures distinctiveness in the reference dataset.

3. **Running YACHT:**
   - This step involves executing the YACHT algorithm to detect the presence of reference genomes within your sample.

### YACHT Commands Overview:

YACHT comprises a suite of commands, each facilitating a specific function in the metagenomic analysis pipeline:

- **yacht_download_demo**: Automatically downloads demo files, enabling you to test the program with pre-set data.
- **yacht_download_genome_sketches**: Retrieves pre-generated Sourmash sketches of reference genomes, ready for use.
- **yacht_sketch_genomes**: Creates Sourmash sketches from your own reference genomes.
- **yacht_training**: Trains models on the sketched reference genomes, preparing them for the detection process.
- **yacht_download_pretrained**: Downloads pre-trained models of reference genomes, ready for immediate use.
- **yacht_sketch_sample**: Generates sketches from your sample's sequencing reads for analysis.
- **yacht_run**: Executes the YACHT algorithm to analyze your sample against the trained reference models.
- **yacht_standardize**: Converts YACHT output into various popular formats for further analysis or reporting.

Each command in YACHT is designed to streamline and simplify the process of metagenomic analysis, from data preparation to final output.

### Creating sketches of your reference database genomes

You will need a reference database in the form of [Sourmash](https://sourmash.readthedocs.io/en/latest/) sketches of a collection of microbial genomes. There are a variety of pre-created databases available at: https://sourmash.readthedocs.io/en/latest/databases.html. Our code uses the "Zipfile collection" format, and we suggest using the [GTDB genomic representatives database](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip).

#### Automatic download of reference sketches

```bash
yacht_download_genome_sketches --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder gtdb-reps
```

#### Manual download of reference sketches

```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
```

#### Automatic sketching of custom database
If you want to use a custom database, you will need to create a Sourmash sketch Zipfile collection from the FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/) for details). In brief, this can be accomplished via the following where your reference fasta files are located in the ref_genomes folder:

```bash
yacht_sketch_genomes --infile ref_genomes --k 31 --scaled 1000 --outfile ref.sig.zip
```

#### Manual sketching of custom database

It also possible to the sketching manully using sourmash directly. 

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


#### Some Pre-trained Reference Databases Available on Zenodo  

For your convenience, we have made available a selection of pre-trained reference databases for GenBank and GTDB genomes on [Zenodo](https://zenodo.org/communities/yacht?q=&l=list&p=1&s=10&sort=newest). These databases can be downloaded and used directly with YACHT, allowing you to skip the training step if they meet your project's requirements.

__Automatic Download__

To automatically download a pre-trained database, use the `yacht_download_pretrained` command as shown below:

```bash
yacht_download_pretrained --database genbank --db_version 2022.3 --ncbi_organism archaea --k 31 --ani_thresh 0.80 --outfolder pretrained_models
```

__Manual Download__

Alternatively, you can manually download a database from Zenodo. Ensure you replace `<zendo_id>` and `<file_name>` with the appropriate values for the database you wish to download:

```bash
# Replace <zendo_id> and <file_name> with your specific values
curl --cookie zenodo-cookies.txt "https://zenodo.org/records/<zendo_id>/files/<file_name>?download=1" --output <file_name>

# Example
curl --cookie zenodo-cookies.txt "https://zenodo.org/records/10113534/files/genbank-2022.03-archaea-k31_0.80_pretrained.zip?download=1" --output genbank-2022.03-archaea-k31_0.80_pretrained.zip
```

### Creating Sketches of Your Sample

Creating a sketch of your sample metagenome is an essential step in the YACHT workflow. This process involves using the same k-mer size and scale factor that were used for the reference database. YACHT provides both automatic and manual methods for sketching.

#### Automatic Sketching of Your Sample

For automated sketching, provide your reads as input. If you have paired-end reads, input both files:

```bash
yacht_sketch_sample --infile <input FASTA/Q file> --k 31 --scaled 1000 --outfile sample.sig.zip
```

#### Manual Method

Alternatively, you can manually sketch your sample using Sourmash commands:

- For a single-end FASTA/Q file:

  ```bash
  sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip <input FASTA/Q file>
  ```

- For paired-end FASTA/Q files, combine them into a single file first:

  ```bash
  cat <FASTA/Q file 1> <FASTA/Q file 2> > combine.fastq (or combine.fasta)
  sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip combine.fastq (or combine.fasta)
  ```

#### Parameters

The Sourmash database offers three k-mer sizes (21, 31, and 51), allowing you to choose the one that best suits your analysis. The scale factor indicates the level of data compression. A smaller scale factor corresponds to a higher proportion of genomes retained in the sketch, which may be preferable for smaller datasets.

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| k=31              | The length of the k-mer.                                     |
| scaled=1000       | The fraction of k-mers to be kept in the sketch.             |
| -o sample.sig.zip | The arbitrary output file name for the sketch files (must be in zip format). |


#### Output

In the two preceding steps, you will obtain a k-mer sketch file in zip format (in each step) for the reference data and the input sample.

| File                      | Content                                                     |
| ------------------------- | ----------------------------------------------------------- |
| gtdb-rs214-reps.k31.zip   | Pre-built k-mer sketch file for GTDB representative genomes |
| training_database.sig.zip | K-mer sketch file for your selected reference genomes       |
| sample.sig.zip            | K-mer sketch file for your input sample                     |




### Preprocess the reference genomes (Training Step)

##### Warning: the training process is time-consuming on large database

In our benchmark with `GTDB representive genomes`, it takes `15 minutes` using `16 threads, 50GB of MEM` on a system equipped with a `3.5GHz AMD EPYC 7763 64-Core Processor`. You can use the pre-trained database (see [here](#some-pre-trained-reference-databases-available-on-zenodo)) to skip this step. The processing time can be significant when executed on GTDB all genomes OR with limited resources. If only part of genomes are needed, one may use `sourmash sig` command to extract signatures of interests only. 


The training script extracts the sketches from the Zipfile-format reference database, and then turns them into a form usable by YACHT. In particular, it removes one of any two organisms that have ANI greater than the user-specified threshold as these two organisms are too close to be "distinguishable".

```bash 
yacht_training--ref_file gtdb-rs214-reps.k31.zip --ksize 31 --num_threads 32 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./
```

#### Parameter

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



### Run the YACHT algorithm

After this, you are ready to perform the hypothesis test for each organism in your reference database. This can be accomplished with something like:

```bash
yacht_run --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --num_threads 32 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result.xlsx
```

#### Parameter

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

### Convert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn)

When we get the EXCEL result file from run_YACHT.py, you can run `standardize_yacht_output.py` under `srcs` folder to covert the YACHT result to other popular output formats (Currently, only `cami`, `biom`, `graphplan` are supported).

__Note__: Before you run `srcs/standardize_yacht_output.py`, you need to prepare a TSV file `genome_to_taxid.tsv` containing two columns: genome ID (genome_id) and its corresponding taxid (taxid). An example can be found [here](demo/toy_genome_to_taxid.tsv). You need to prepare it according to the reference database genomes you used. 

Then you are ready to run `standardize_yacht_output.py` with something like:
```bash
yacht_standardize --yacht_output 'result.xlsx' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir ./
```

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --yacht_output    | the path to the output excel file generated by `run_YACHT.py` |
| --sheet_name      | specify which spreadsheet result you want to covert from     |
| --genome_to_taxid | the path to the location of `genome_to_taxid.tsv` you prepared |
| --mode            | specify to which output format you want to convert (e.g., 'cami', 'biom', 'graphplan')
| --sample_name     | A random name you would like to show in header of the cami file. Default: Sample1.' |
| --outfile_prefix  | the prefix of the output file. Default: result | 
| --outdir          | the path to output directory where the results will be genreated |


