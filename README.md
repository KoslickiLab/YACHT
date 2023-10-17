# YACHT

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity.

The associated preprint can be found at:  https://doi.org/10.1101/2023.04.18.537298. Please cite via:

>Koslicki, D., White, S., Ma, C., & Novikov, A. (2023). YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample. bioRxiv, 2023-04.

</br>

## Quick start

```
cd demo

# build k-mer sketch for the query data and ref genomes
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip  query_data/query_data.fq
sourmash sketch fromfile ref_paths.csv -p dna,k=31,scaled=1000,abund -o ref.sig.zip 

# prepare reference k-mer dictionary its k-mer sketch
python ../make_training_data_from_sketches.py --ref_file ref.sig.zip --ksize 31 --out_prefix 'demo_ani_thresh_0.95' --ani_thresh 0.95

# run YACHT algorithm to check the presence of reference genomes in the input sample
python ../run_YACHT.py --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --min_coverage 1 0.6 0.2 0.1 --outdir './'

# convert result to CAMI profile format (TBD)
python ../srcs/standardize_yacht_output.py --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name MySample --outfile_prefix cami_result --outdir './'
```

There will be an output EXCEL file `result.xlsx` recoding the presence of reference genomes with the given minimum coverage of `1 0.6 0.2 0.1`

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

YACHT requires Python 3 or higher. We recommend using a virtual environment (such as conda) to install YACHT. To create a virtual environment, run:

```bash
# Clone the repo
git clone https://github.com/KoslickiLab/YACHT.git
cd YACHT

# Create a conda environment for YACHT
conda env create -f env/yacht_env.yaml
conda activate yacht
```

</br>

## Usage

The workflow for YACHT is as follows: 

1. Create sketches of your reference database genomes and of your sample
2. Create a reference dictionary matrix of the reference sketch
3. Run YACHT to detect the presence of reference genomes in your sample

</br>

### Creating sketches of your reference database genomes

You will need a reference database in the form of [Sourmash](https://sourmash.readthedocs.io/en/latest/) sketches of a collection of microbial genomes. There are a variety of pre-created databases available at: https://sourmash.readthedocs.io/en/latest/databases.html. Our code uses the "Zipfile collection" format, and we suggest using the [GTDB genomic representatives database](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip):

```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.zip
```

If you want to use a custom database, you will need to create a Sourmash sketch of your FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/) for details). In brief, this can be accomplished via the following:

If you have a single FASTA file with _one genome_ per record:

```bash
sourmash sketch dna -f -p k=31,scaled=1000,abund --singleton <your multi-FASTA file> -o training_database.sig.zip
```

If you have a directory of FASTA files, one per genome:

```bash
# cd into the relevant directory
sourmash sketch dna -f -p k=31,scaled=1000,abund *.fasta -o ../training_database.sig.zip
# cd back to YACHT
```

</br>

### Creating sketches of your sample

You will then create a sketch of your sample metagenome, using the same k-mer size and scale factor

```bash
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip <input FASTA/Q file>
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

### Creating a reference dictionary matrix

The script `make_training_data_from_sketches.py` collects and transforms the sketched microbial genomes, getting them into a form usable by YACHT. In particular, it removes one of any two organisms that are within the ANI threshold the user specifies as making two organisms "indistinguishable".

```bash 
python make_training_data_from_sketches.py --ref_file 'gtdb-rs214-reps.k31.zip' --ksize 31 --out_prefix 'gtdb_ani_thresh_0.95' --ani_thresh 0.95
```

#### Parameter

The most important parameter of this command is `--ani_thresh`: this is average nucleotide identity (ANI) value below which two organisms are considered distinct. For example, if `--ani_thresh` is set to 0.95, then two organisms with ANI >= 0.95 will be considered indistinguishable. Only the largest of such organisms will be kept in the reference dictionary matrix. The default value of `--ani_thresh` is 0.95. The `--ani_thresh` value chosen here must match the one chosen for the YACHT algorithm (see below).  

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --ani_thresh 0.95 | the cutoff by which two organisms are considered indistinguishable |
| --ksize 31        | the length of k-mer, must match the k size used in previous sketching steps |
| --out_prefix      | prefix for output files (see details below)                  |

#### Output (to check after Chunyu's update)

| File (names starting with out_prefix) | Content                                                      |
| ------------------------------------- | ------------------------------------------------------------ |
| _config.json                          | A json file stores input files needed to run the next YACHT algorithm |
| _hash_to_col_idx                      | A pickle object storing the the reference dictionary matrix  |
| _processed_org_idx.csv                | A csv file storing statistics of reference genomes           |
| _ref_matrix_processed.npz             |                                                              |



</br>

### Run the YACHT algorithm

After this, you are ready to perform the hypothesis test for each organism in your reference database. This can be accomplished with something like:

```bash
python run_YACHT.py --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --outdir './'
```

#### Parameter

The `--significance` parameter is basically akin to your confidence level: how sure do you want to be that the organism is present? Higher leads to more false negatives, lower leads to more false positives. 
The `--min_coverage` parameter dictates what percentage (value in `[0,1]`) of the distinct k-mers (think: whole genome) must have been sequenced and present in my sample to qualify as that organism as being "present." Setting this to 1 is usually safe, but if you have a very low coverage sample, you may want to lower this value. Setting it higher will lead to more false negatives, setting it lower will lead to more false positives (pretty rapidly).

| Parameter                               | Explanation                                                  |
| --------------------------------------- | ------------------------------------------------------------ |
| --json gtdb_ani_thresh_0.95_config.json | This file contains results from the previous step generating reference dictionary matrix |
| --significance                          | Statistical confidence level                                 |
| --min_coverage                          | minimum genome coverage for an organism to be "present"      |

#### Output

The output file will be a CSV file; column descriptions can be found [here](docs/column_descriptions.csv). The most important are the following:

* `organism`: The name of the organism
* `in_sample_est`: This value is either 0 or 1: if 0, there was not enough evidence to claim this organism is present in the sample. 
* `p_vals`: Probability of observing this or more extreme result at the given ANI threshold, assuming the null hypothesis.

Other interesting columns include:

* `num_exclusive_kmers`: How many k-mers were found in this organism and no others
* `num_matches`: How many k-mers were found in this organism and the sample
* `acceptance_threshold_*`: How many k-mers must be found in this organism to be considered "present" at the given ANI threshold. Hence, `in_sample_est` is 1 if `num_matches` >= `acceptance_threshold_*` (adjusting by coverage if desired).
* `alt_confidence_mut_rate`: What the mutation rate (1-ANI) would need to be to get your false positive to match the false negative rate of 1-`significance`.

</br>

### Convert YACHT result to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPlAn)

When we get the EXCEL result file from run_YACHT.py, you run `standardize_yacht_output.py` to covert the YACHT result to other popular output formats. Currently, only `cami`, `biom`, `graphplan` are supported. (__Note: Before you run `srcs/standardize_yacht_output.py`, you need to first prepare a file `genome_to_taxid.tsv` which is a TSV file with two columns: genome ID (genome_id) and its corresponding taxid (taxid)__).

```bash
python srcs/standardize_yacht_output.py --yacht_output 'result.xlsx' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir './'
```

Note: we may need to build a GTDB-taxid map.

```
# error message:
pytaxonkit.TaxonKitCLIError: 15:11:31.257 [ERRO] taxonomy data not found, please download and uncompress ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz, and copy "names.dmp", "nodes.dmp", "delnodes.dmp", and "merged.dmp" to /home/grads/sml6467/.taxonkit
```


