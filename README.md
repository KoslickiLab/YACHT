# YACHT

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).

The associated preprint can be found at:  https://doi.org/10.1101/2023.04.18.537298. Please cite via:

>Koslicki, D., White, S., Ma, C., & Novikov, A. (2023). YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. bioRxiv, 2023-04.

</br>

## Quick start
We provide a demo to show how to use YACHT. Please follow the command lines below to try it out:

```bash
NUM_THREADS=64 # if your machine doesn't have so many CPU cores, feel free to reduce this value.

cd demo

# build k-mer sketches for the query sample and ref genomes
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip query_data/query_data.fq
sourmash sketch fromfile ref_paths.csv -p dna,k=31,scaled=1000,abund -o ref.sig.zip 

# preprocess the reference genomes (training step)
python ../make_training_data_from_sketches.py --ref_file ref.sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./

# run YACHT algorithm to check the presence of reference genomes in the query sample (inference step)
python ../run_YACHT.py --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads ${NUM_THREADS} --min_coverage_list 1 0.6 0.2 0.1 --outdir ./

# convert result to CAMI profile format (Optional)
python ../srcs/standardize_yacht_output.py --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./
```

There will be an output EXCEL file `result.xlsx` recoding the presence of reference genomes with different spreadsheets given the minimum coverage of `1 0.6 0.2 0.1`.

</br>

### Contents

- [Installation](#installation)
  * [Conda](#conda)
  * [Manual installation](#manual-installation)
- [Usage](#usage)
  * [Creating sketches of your reference database genomes](#creating-sketches-of-your-reference-database-genomes)
  * [Creating sketches of your sample](#creating-sketches-of-your-sample)
    + [Parameters](#parameters)
    + [Output](#output)
  * [Creating a reference dictionary matrix](#creating-a-reference-dictionary-matrix)
    + [Parameter](#parameter)
    + [Output (to check after Chunyu's update)](#output-to-check-after-chunyus-update)
  * [Run the YACHT algorithm](#run-the-yacht-algorithm)
    + [Parameter](#parameter-1)
    + [Output](#output-1)
  * [Convert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn)](#convert-yacht-result-to-other-popular-output-formats-eg-cami-profiling-format-biom-format-graphplan)


## Installation

### Conda

A conda release will be coming soon. In the meantime, please install manually.

### Manual installation

YACHT requires Python 3 or higher. We recommend using a virtual environment (such as [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)) to run YACHT. To create a virtual environment, run:

```bash
# Clone the repo
git clone https://github.com/KoslickiLab/YACHT.git
cd YACHT

# Set up an environment for YACHT
bash setup.sh

# Activiate YACHT environment
conda activate yacht_env
```

</br>

## Usage

The workflow for YACHT is as follows: 

1. Create sketches of your reference database genomes and of your sample
2. Preprocess the reference genomes by removing the "too similiar" genomes based on `ANI` using the `ani_thresh` parameter 
3. Run YACHT to detect the presence of reference genomes in your sample

### Creating sketches of your reference database genomes

You will need a reference database in the form of [Sourmash](https://sourmash.readthedocs.io/en/latest/) sketches of a collection of microbial genomes. There are a variety of pre-created databases available at: https://sourmash.readthedocs.io/en/latest/databases.html. Our code uses the "Zipfile collection" format, and we suggest using the [GTDB genomic representatives database](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip):

```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
```

If you want to use a custom database, you will need to create a Sourmash sketch Zipfile collection from the FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/) for details). In brief, this can be accomplished via the following:

If you have a single FASTA file with _one genome_ per record:

```bash
sourmash sketch dna -f -p k=31,scaled=1000,abund --singleton <your multi-FASTA file> -o training_database.sig.zip
```

If you have a directory of FASTA files, one per genome:

```bash
## Method 1
# cd into the relevant directory
sourmash sketch dna -f -p k=31,scaled=1000,abund *.fasta -o ../training_database.sig.zip
# cd back to YACHT

## Method 2 
# put all full paths of FASTA/FASTQ file into a file, one path per line
find <path of foler containg FASTA/FASTQ files> > dataset.csv
sourmash sketch fromfile dataset.csv -p dna,k=31,scaled=1000,abund -o ../training_database.sig.zip
# cd back to YACHT
```

### Creating sketches of your sample

You will then create a sketch of your sample metagenome, using the same k-mer size and scale factor

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

### Preprocess the reference genomes (Training Step)

The script `make_training_data_from_sketches.py` extracts the sketches from the Zipfile-format reference database, and then turns them into a form usable by YACHT. In particular, it removes one of any two organisms that have ANI greater than the user-specified threshold as these two organisms are too close to be "distinguishable".

```bash 
python make_training_data_from_sketches.py --ref_file gtdb-rs214-reps.k31.zip --ksize 31 --num_threads 32 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./
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
| _manifest.tsv                         | A TSV file contains organisms and their relevant info after removing the similiar ones.
| _rep_to_corr_orgas_mapping.tsv       | A TSV file contains representaive organisms and their similiar organisms that have been removed |


</br>

### Run the YACHT algorithm

After this, you are ready to perform the hypothesis test for each organism in your reference database. This can be accomplished with something like:

```bash
python run_YACHT.py --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --num_threads 32 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --outdir ./
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
| --out_filename                          | filename of output excel result (default: 'result.xlsx') |
| --outdir                          | the path to output directory where the results and intermediate files will be genreated |

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
python srcs/standardize_yacht_output.py --yacht_output 'result.xlsx' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir ./
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



