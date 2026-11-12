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
yacht run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --min_coverage_list 1 0.6 0.2 0.1 --outdir ./

# convert result to CAMI profile format (Optional)
yacht convert --yacht_output_dir ./results --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./
```

The output will be stored in the `results` folder containing:
- `result.xlsx`: An EXCEL file recording the presence of reference genomes with different spreadsheets given the minimum coverage of `1 0.6 0.2 0.1`.
- `result_all.txt`: A TXT file containing all unfiltered results for all user-given min_coverage values.

> **New in YACHT 2.0:** `yacht train` also writes a binary index of the reference database (`<prefix>_index.yidx`) that makes `yacht run` much faster (it no longer parses every reference genome sketch). `yacht run` is single-threaded (its `--num_threads` option is deprecated and ignored), so many samples can be processed in parallel. Databases created with older YACHT versions need the index to be added once, see [here](#using-databases-created-with-yacht--20).

</br>

### Contents

- [Quick demonstration](#quick-demonstration)
- [Installation](#installation)
  - [Conda installation](#conda-installation)
  - [Manual installation](#manual-installation)
    - [Using Conda](#using-conda)
    - [Using Mamba](#using-mamba)
    - [Using Docker](#using-docker)
- [Commands](#commands)
- [Pretrained reference databases](#pretrained-reference-databases)
- [Workflow](#workflow)
  - [1. Create sketches of reference database genomes and samples](#1-create-sketches-of-reference-database-genomes-and-samples)
    - [Reference sketches](#reference-sketches)
    - [Sample sketches](#sample-sketches)
  - [2. Preprocess the reference genomes (yacht train)](#2-preprocess-the-reference-genomes-yacht-train)
    - [Parameters](#parameters)
    - [Output](#output)
    - [Using databases created with YACHT < 2.0](#using-databases-created-with-yacht--20)
  - [3. Run the YACHT algorithm (yacht run)](#3-run-the-yacht-algorithm-yacht-run)
    - [Parameters](#parameters-1)
    - [Output](#output-1)
  - [4. Convert YACHT result to other popular output formats (yacht convert)](#4-convert-yacht-result-to-other-popular-output-formats-yacht-convert)
    - [Parameters](#parameters-2)

## Installation

### Conda installation

YACHT is available on Bioconda and can be installed with:
```bash
# create conda environment
conda create -n yacht_env

# activate environment
conda activate yacht_env

# install YACHT
conda install -c conda-forge -c bioconda yacht
```

### Manual installation
YACHT requires **Python >3.6** (and <3.12) with the following core genomics dependencies: `sourmash` (>=4.8.3), `sourmash_plugin_branchwater`, and `pytaxonkit`. Installing from source also requires a C++ compiler (with `make`) and a [Rust](https://www.rust-lang.org/tools/install) toolchain (`cargo`, >=1.80), which are used to compile the native parts of YACHT (the `yacht train` comparison core and the `yacht run` index engine); the conda environment file below installs Rust for you. The full list of dependencies can be found in the [environment configuration](https://github.com/KoslickiLab/YACHT/blob/main/env/yacht_env.yml). To ensure a clean and isolated workspace, we recommend using a virtual environment. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba), a faster alternative to Conda.
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
If you prefer using Mamba instead of Conda, simply replace `conda` with `mamba` in the above commands.

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
  # Example for downloading the k31 sketches of representative genomes of GTDB rs226 version 
  yacht download default_ref_db --database gtdb --db_version rs226 --gtdb_type reps --k 31 --outfolder ./
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | database          | two options for default reference databases: 'genbank' or 'gtdb' |
  | db_version        | the version of database, options: "genbank-2022.03", "genbank-2025.04" (fungi and protozoa only), "rs202", "rs207", "rs214", "rs220", "rs226" |
  | ncbi_organism     | the NCBI organism for the NCBI reference genome, options: "archaea", "bacteria", "fungi", "virus", "protozoa"|
  | gtdb_type         | for GTDB database, chooses "reps" (representative genomes) or "full" (all genomes; the default) |
  | k                 | the length of k-mer |
  | outfolder         | the path to a folder where the downloaded file is expected to locate |


  + `pretrained_ref_db` can automatically download our pre-trained reference genome databases (see [Pretrained reference databases](#pretrained-reference-databases)), which can be used directly as input for the `yacht run` module. The available databases are looked up on [Zenodo](https://zenodo.org/communities/yacht) when the command runs, so databases added there later can be downloaded without updating YACHT.
  ```bash
  # List the available pretrained databases (name, ANI thresholds, size)
  yacht download pretrained_ref_db --list
  # Example: the database trained from the GTDB R232 representative genomes with k=31 and ani_threshold=0.95
  yacht download pretrained_ref_db --name gtdb-rs232-reps.k31 --ani_thresh 0.95 --outfolder ./
  # The same database selected by its family and version
  yacht download pretrained_ref_db --database gtdb --db_version rs232 --k 31 --ani_thresh 0.95 --outfolder ./
  ```
  | Parameter         | Explanation                                                  |
  | ----------------- | ------------------------------------------------------------ |
  | list              | list the pretrained databases available on Zenodo and exit |
  | name              | the name of a pretrained database, as shown by `--list` (instead of `--database`, `--db_version`, `--ncbi_organism` and `--k`) |
  | database          | the database family, e.g. 'genbank' or 'gtdb' |
  | db_version        | the version of database, e.g. "genbank-2022.03" or "rs214" |
  | ncbi_organism     | the NCBI organism for the NCBI reference genome, options: "archaea", "bacteria", "fungi", "virus", "protozoa"|
  | ani_thresh      | the ANI threshold the database was trained with (see `--list`) |
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

- The `train` module preprocesses the given sketches of reference genomes (the `.zip` file) to identify and merge the "identical" genomes based on the given ANI threshold (e.g., `--ani_thresh 0.95`). For an example, please refer to the `yacht train` command in the [Quick demonstration](#quick-demonstration) section.

- The `run` module runs the YACHT algorithm to detect the presence of reference genomes in a given sample. For an example, please refer to the `yacht run` command in the [Quick demonstration](#quick-demonstration) section.

- The `convert` module can convert YACHT results to other popular output formats (e.g., CAMI profiling format, BIOM format, GraphPhlAn). For an example, please refer to the `yacht convert` command in the [Quick demonstration](#quick-demonstration) section.

## Pretrained reference databases

Pretrained databases can be used with `yacht run` directly, without the training step. They are published in the [YACHT community on Zenodo](https://zenodo.org/communities/yacht), where `yacht download pretrained_ref_db` finds them:

```bash
yacht download pretrained_ref_db --list
yacht download pretrained_ref_db --name gtdb-rs232-reps.k31 --ani_thresh 0.95 --outfolder ./
yacht run --json gtdb-rs232-reps.k31_0.95_pretrained/gtdb-rs232-reps.k31_0.95_config.json --sample_file sample.sig.zip --outdir ./
```

The databases built for YACHT 2.0 are trained at k=31 with the ANI thresholds 0.80, 0.95, 0.995 and 0.9995. Each one comes with a data card (the `README.md` in the zip and on its Zenodo record): sources, licenses, citations, identifiers, genome counts, and how coarse the database is at its ANI threshold.

| database | content | genomes (before training) | scaled | data cards (ANI) |
|---|---|---:|---:|---|
| `gtdb-rs232-reps.k31` | GTDB R232 species representatives (bacteria, archaea) | 199,923 | 1000 | [0.80](docs/database_cards/gtdb-rs232-reps.k31_0.80.md) · [0.95](docs/database_cards/gtdb-rs232-reps.k31_0.95.md) · [0.995](docs/database_cards/gtdb-rs232-reps.k31_0.995.md) · [0.9995](docs/database_cards/gtdb-rs232-reps.k31_0.9995.md) |
| `globdb-r232.k31` | GlobDB r232: GTDB R232 species representatives + species-level MAGs of 25 catalogues | 346,233 | 1000 | [0.80](docs/database_cards/globdb-r232.k31_0.80.md) · [0.95](docs/database_cards/globdb-r232.k31_0.95.md) · [0.995](docs/database_cards/globdb-r232.k31_0.995.md) · [0.9995](docs/database_cards/globdb-r232.k31_0.9995.md) |
| `genbank-2025.04-fungi-k31` | GenBank fungi, April 2025, without the assemblies NCBI flags as contaminated | 21,070 | 1000 | [0.80](docs/database_cards/genbank-2025.04-fungi-k31_0.80.md) · [0.95](docs/database_cards/genbank-2025.04-fungi-k31_0.95.md) · [0.995](docs/database_cards/genbank-2025.04-fungi-k31_0.995.md) · [0.9995](docs/database_cards/genbank-2025.04-fungi-k31_0.9995.md) |
| `genbank-2025.04-protozoa-k31` | GenBank protozoa, April 2025, without the assemblies NCBI flags as contaminated | 2,346 | 1000 | [0.80](docs/database_cards/genbank-2025.04-protozoa-k31_0.80.md) · [0.95](docs/database_cards/genbank-2025.04-protozoa-k31_0.95.md) · [0.995](docs/database_cards/genbank-2025.04-protozoa-k31_0.995.md) · [0.9995](docs/database_cards/genbank-2025.04-protozoa-k31_0.9995.md) |
| `ictv-vmr-msl40.k31-sc50` | ICTV exemplar viruses (VMR MSL40) | 17,765 | 50 | [0.80](docs/database_cards/ictv-vmr-msl40.k31-sc50_0.80.md) · [0.95](docs/database_cards/ictv-vmr-msl40.k31-sc50_0.95.md) · [0.995](docs/database_cards/ictv-vmr-msl40.k31-sc50_0.995.md) · [0.9995](docs/database_cards/ictv-vmr-msl40.k31-sc50_0.9995.md) |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | NCBI RefSeq viral genomes, January 2025 | 16,280 | 50 | [0.80](docs/database_cards/ncbi-viruses-refseq-2025.01.k31-sc50_0.80.md) · [0.95](docs/database_cards/ncbi-viruses-refseq-2025.01.k31-sc50_0.95.md) · [0.995](docs/database_cards/ncbi-viruses-refseq-2025.01.k31-sc50_0.995.md) · [0.9995](docs/database_cards/ncbi-viruses-refseq-2025.01.k31-sc50_0.9995.md) |
| `ncbi-viruses-2025.01.k31-sc50` | all NCBI viral assemblies, January 2025 (also GenBank assemblies that are not single viral genomes, see its data card) | 230,973 | 50 | [0.80](docs/database_cards/ncbi-viruses-2025.01.k31-sc50_0.80.md) · [0.95](docs/database_cards/ncbi-viruses-2025.01.k31-sc50_0.95.md) · [0.995](docs/database_cards/ncbi-viruses-2025.01.k31-sc50_0.995.md) · [0.9995](docs/database_cards/ncbi-viruses-2025.01.k31-sc50_0.9995.md) |
| `yacht-union-insdc-2026.k31` | bacteria, archaea, fungi, protozoa and viruses (GTDB R232 + fungi + protozoa + RefSeq and ICTV viruses), NCBI accessions only: easiest to go from a result to the genome | 257,384 | 1000 | [0.80](docs/database_cards/yacht-union-insdc-2026.k31_0.80.md) · [0.95](docs/database_cards/yacht-union-insdc-2026.k31_0.95.md) · [0.995](docs/database_cards/yacht-union-insdc-2026.k31_0.995.md) · [0.9995](docs/database_cards/yacht-union-insdc-2026.k31_0.9995.md) |
| `yacht-union-globdb-2026.k31` | the same with GlobDB r232 instead of GTDB R232: most comprehensive | 403,694 | 1000 | [0.80](docs/database_cards/yacht-union-globdb-2026.k31_0.80.md) · [0.95](docs/database_cards/yacht-union-globdb-2026.k31_0.95.md) · [0.995](docs/database_cards/yacht-union-globdb-2026.k31_0.995.md) · [0.9995](docs/database_cards/yacht-union-globdb-2026.k31_0.9995.md) |

[docs/databases.md](docs/databases.md) explains which database to choose, the genome identifiers and taxonomy, known limitations and how the databases were built; [docs/retrieving_genomes.md](docs/retrieving_genomes.md) explains how to get the genome sequence of a result. All data cards are listed in [docs/database_cards](docs/database_cards/README.md).

- **YACHT version.** The databases built for YACHT 2.0 need YACHT 2.0.0 or later: they hold the YACHT 2.0 index of the genomes but not the genome sketches that YACHT 1.x reads. Their config file records this (`min_yacht_version`, and the YACHT version that made them as `yacht_version`), so that YACHT refuses a database that needs a newer version, and their data cards and Zenodo records say so. The older databases in the Zenodo community (GenBank 2022.03, GTDB rs214), made with YACHT 1.x, also work with YACHT 2.x: `yacht download pretrained_ref_db` adds the index (see [Using databases created with YACHT < 2.0](#using-databases-created-with-yacht--20)).
- **ANI 0.80 is very coarse,** much coarser than genus level (at k=31 it amounts to a k-mer containment of 0.001): in `gtdb-rs232-reps.k31_0.80`, 70% of the GTDB genera keep no genome of their own. Use 0.95 or higher for genus- or species-level results; each data card has a table of how coarse its database is.
- **Viruses** are small, so the viral databases use scaled 50: sketch the sample with `yacht sketch sample --scaled 50`. A sample sketched with a smaller scaled than a database is downsampled to it, so the same sketch also works with the scaled-1000 databases.
- **Licenses:** CC BY-SA 4.0 for the databases built from GTDB or GlobDB (including the unions), CC BY 4.0 for the ICTV VMR, CC0 for the NCBI-based ones; the data card of each database lists what to cite.
- To download a database without YACHT, get `<name>_<ANI>_pretrained.zip` from its Zenodo record and unzip it; `yacht run` finds the database files next to the config file.

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
#### Reference sketches


Use the command `yacht sketch` to generate sketches for both the samples and the reference genomes. Users must utilize sourmash to extract sketches from a reference database of microbial genomes. [sourmash Databases](https://sourmash.readthedocs.io/en/latest/databases.html) provide a variety of pre-formed databases of such sketches, or users can create a custom database using the `sourmash sketch` command on FASTA/FASTQ files of reference genomes (see the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/)). Other available databases include the [GTDB genomic representatives database](https://gtdb.ecogenomic.org/downloads). The sketches for samples must be generated using the same $k$-mer size and scale factor as those used for the reference database. The scale factor acts as an indicator of data compression, with smaller values being more appropriate for smaller datasets. 

We suggest trying the reference sketches prepared by the sourmash team, e.g. of the [GTDB rs226 representative genomes](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226-reps.k31.sig.zip):

```bash
yacht download default_ref_db --database gtdb --db_version rs226 --gtdb_type reps --k 31 --outfolder ./
```

Or
```bash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/gtdb-rs226-reps.k31.sig.zip
```

For custom databases, you will need to create a Sourmash sketch Zipfile collection from the FASTA/FASTQ files of your reference database genomes (see [Sourmash documentation](https://sourmash.readthedocs.io/en/latest/)). Following commands accomplish it:

A single FASTA file with _one genome_ per record:

```bash
# This is equivalent to: sourmash sketch dna -f -p k=31,scaled=1000,abund --singleton <path to your multi-FASTA file> -o training_database.sig.zip
yacht sketch ref --infile <path to your multi-FASTA file> --kmer 31 --scaled 1000 --outfile training_database.sig.zip
```

A directory of FASTA files, one per genome:

```bash
# This is equivalent to: find <path of folder containing FASTA/FASTQ files> > dataset.csv; sourmash sketch fromfile dataset.csv -p dna,k=31,scaled=1000,abund -o training_database.sig.zip
yacht sketch ref --infile <path of folder containing FASTA/FASTQ files> --kmer 31 --scaled 1000 --outfile training_database.sig.zip
```

</br>


#### Sample sketches

This process should use the same k-mer size and scale factor that were used for the reference database. 

```bash
# For a single-end FASTA/Q file (plain or gzipped)
# the command below gives the same sketch (hashes and abundances) as: sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip <input FASTA/Q file>
yacht sketch sample --infile <input FASTA/Q file> --kmer 31 --scaled 1000 --outfile sample.sig.zip

# For pair-end FASTA/Q files, you need to separately specify two FASTA/Q files
# the command below gives the same sketch (hashes and abundances) as: cat <FASTA/Q file 1> <FASTA/Q file 2> > combine.fastq (or combine.fasta); sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip combine.fastq (or combine.fasta)
yacht sketch sample --infile <FASTA/Q file 1> <FASTA/Q file 2> --kmer 31 --scaled 1000 --outfile sample.sig.zip
```

> **New in YACHT 2.0:** `yacht sketch sample` uses `sourmash scripts singlesketch` from the sourmash branchwater plugin (`sourmash_plugin_branchwater` 0.9.13 or later). It computes the sketch with sourmash's own code on one CPU core, without per-read Python overhead and without first copying paired-end files into one file. On the SRA metagenomes used for the [`yacht run` timings](#3-run-the-yacht-algorithm-yacht-run) (2.2–3.1 Gbp of paired-end FASTQ), sketching took 38–60 s instead of 88–140 s, and 0.3 GB of memory instead of up to 9.6 GB. Gzipped paired-end files are now read directly. The sketch is named after the first input file (before, it had no name). With an older plugin, or an output file that does not end in `.zip`, `.sig` or `.sig.gz`, `yacht sketch sample` uses `sourmash sketch dna` as before. Sketching is still slower than `sylph sketch` (4–6 s at c1000 and 8–13 s at c200 on the same samples), so from raw reads to a profile, sylph is faster.

Note: Sourmash database offers three available k values (21, 31, and 51), allowing you to select the one that best suits your particular analytical needs. The scale factor serves as an indicator of data compression, and if your dataset is small, you might consider using a smaller value (corresponding to a higher portion of genomes retained in the sketch).


</br>

### 2. Preprocess the reference genomes (yacht train)

`yacht train` identifies and merges genomes that are roughly identical based on Average Nucleotide Identity (ANI). The module utilizes a fast algorithm written by C++ to preprocess the reference genomes. It then builds a binary index of the selected genomes (`<prefix>_index.yidx`) that `yacht run` uses (building it for the ~80,000 GTDB r214 representative genomes takes under a minute and ~3 GB of RAM). In our test with the GTDB representative genomes (r214) including `85,205` species-level genomes, YACHT takes around `12 minutes` and `52 GB` of RAM to preprocess them and generate the reference files on a Ubuntu 22.04.5 system using 64 threads; the `199,923` GTDB R232 representative genomes took 18–22 minutes and 84 GB with 48 threads (the costs of all pretrained databases are in [databases/README.md](databases/README.md#cost-measured-on-a-768-core--3-tb-server-october-2026)). You can also use the [pretrained databases](#pretrained-reference-databases) we built to skip this step.


```bash 
yacht train --ref_file gtdb-rs226-reps.k31.sig.zip --ksize 31 --num_threads 64 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./
```

#### Parameters

The most important parameter of this script is `--ani_thresh`: this is average nucleotide identity (ANI) value equal to or below which two organisms are considered distinct. For example, if `--ani_thresh` is set to 0.95, then two organisms with ANI > 0.95 will be considered indistinguishable. For the organisms with ANI > 0.95, only the one with the largest number of unique kmers will be kept. If there is a tie in the number of unique kmers, one organism will be randomly selected. The default value of `--ani_thresh` is 0.95. The ANI threshold is converted to a k-mer containment threshold (ANI^k) for the comparison of the genome sketches; low thresholds are coarse (at k=31, 0.80 is much coarser than genus level, see [Pretrained reference databases](#pretrained-reference-databases)).  

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --ref_file        | the path to the sourmash signature database zip file |
| --ksize           | the length of k-mer, must match the k size used in previous sketching steps (default: 31) |
| --num_threads     | the number of threads to use for parallelization (default: 16) |
| --ani_thresh      | the cutoff by which two organisms are considered indistinguishable, greater than 0 and at most 1 (default: 0.95) |
| --prefix          | the prefix for output folders and files (see details below) |
| --outdir          | the path to output directory where the results and intermediate files will be generated |

#### Output

| File (names starting with prefix)     | Content                                                      |
| ------------------------------------- | ------------------------------------------------------------ |
| _config.json                          | A JSON file storing the information `yacht run` needs (paths, k-mer size, scaled, ANI threshold) and the YACHT version that made the database (`yacht_version`) |
| _processed_manifest.tsv               | A TSV file contains organisms and their relevant info after removing the similar ones |
| _index.yidx                           | A binary index of the k-mers of the selected organisms, used by `yacht run` (new in YACHT 2.0) |
| _intermediate_files/                  | A folder with the signatures of the reference genomes and other intermediate files |

> **New in YACHT 2.0:** the genomes are identified by their signature names in all outputs, so `yacht train` stops with an error if two reference sketches have the same name (sketches made with `sourmash sketch dna` without `--name` have no name; `yacht sketch ref` names them by their file names). Empty reference sketches (genomes too short for any hash at this k-mer size and scaled) are skipped with a warning, since they can never be detected, a zip holding sketches of several k-mer sizes is trained on the `--ksize` sketches only, and zips holding identical sketches (which sourmash stores as `signatures/<md5>.sig.gz_0`, ...) are read correctly; YACHT 1.3.x crashed in these four cases. Sketches that are not valid FracMinHash sketches are rejected (see [below](#3-run-the-yacht-algorithm-yacht-run)). For large databases, the genome comparison now runs in several passes, so that its comparison matrix stays within ~32 GB of memory (it used to need 4 bytes per pair of genomes, e.g. 160 GB for 200,000 genomes); this does not change the result.

#### Using databases created with YACHT < 2.0

Since YACHT 2.0, `yacht run` requires the binary index `<prefix>_index.yidx` next to `<prefix>_config.json`. Databases trained with older YACHT versions (and pretrained databases downloaded manually or with older YACHT versions) do not have it, and `yacht run` stops with instructions. The index can be added in place, without the original reference sketches, with:

```bash
python -m yacht.yacht_index --json gtdb_ani_thresh_0.95_config.json
```

This reads the genome signatures in `<prefix>_intermediate_files/signatures`, writes `<prefix>_index.yidx` and records its location in the config file (as `index_file_path`; older YACHT versions ignore this entry, so the database stays usable with them). Building the index for ~80,000 GTDB genomes takes under a minute. The index is tied to the manifest it was built for: if the manifest changes, `yacht run` reports the index as out of date, and `python -m yacht.yacht_index --json <config> --force` rebuilds it. To distribute a pretrained database that includes the index, simply zip the database folder after running the command above (`yacht download pretrained_ref_db` uses an index found in the downloaded folder as is).

The other way round does not work: YACHT 1.x cannot use the [pretrained databases](#pretrained-reference-databases) built for YACHT 2.0, which leave out the genome sketches (their config records `min_yacht_version`). A database trained with `yacht train` 2.x keeps its genome sketches in `<prefix>_intermediate_files` and remains usable with YACHT 1.x.


</br>

### 3. Run the YACHT algorithm (yacht run)

After this, you are ready to perform the hypothesis test via `yacht run` for each organism in your reference database. This can be accomplished with something like:

```bash
yacht run --json 'gtdb_ani_thresh_0.95_config.json' --sample_file 'sample.sig.zip' --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --outdir ./
```

`yacht run` memory-maps the database index (`<prefix>_index.yidx`) instead of reading the genome signatures, and runs on a single CPU core. With the pretrained GTDB rs214 representative databases (78,861 genomes at ANI 0.95 or 85,205 at ANI 0.9995; ~2 GB index), a sketched short-read metagenome (human gut, soil and marine SRA runs) took 1.1–2.5 s on one CPU core, including Python start-up and reading the index from disk. `sylph profile -t 1` took 2.5–5.5 s on the same samples with Sylph's c1000 GTDB r214 database and 15–22 s with its default c200 database. Sketching the sample is not included in these times (see [sample sketches](#sample-sketches)). Peak memory is ~2 GB, almost all of it the memory-mapped index. With the larger YACHT 2.0 pretrained databases, e.g. GTDB R232 at ANI 0.95 (184,616 genomes), the same samples took 2–3 s and ~5 GB. Because the index is memory-mapped, concurrent `yacht run` processes (e.g. one per sample) share a single copy of it in memory. To process many samples, run one `yacht run` per sample in parallel, using a different `--outdir` for each sample.

> **New in YACHT 2.0:**
> - A sample sketched with a smaller `--scaled` than the reference database (e.g. 100 for a database built with 1000) is downsampled to the database's scaled, which gives exactly the sketch the sample would have at that scaled. One sample sketch can thus be used with databases of different scaled values. A sample sketched with a larger scaled than the database is an error.
> - A FracMinHash sketch is a set of distinct hash values, each at most the sketch's `max_hash`, with an abundance of at least 1 for each hash. sourmash only writes such sketches; a sample sketch that is not one (repeated hash values, hashes above its `max_hash`, abundances of 0 or not one per hash) is rejected with an error.
> - The column `num_exclusive_kmers_in_sample_sketch` holds the number of distinct k-mers (hashes) in the sample sketch, as [documented](docs/column_descriptions.csv). YACHT 1.3.x reported the sample's mean k-mer abundance in this column (it was not used in any computation).
> - If a sheet of `result.xlsx` would have more rows than Excel allows (1,048,575 organisms), `result.xlsx` is not written and `yacht run` stops with an error; `result_all.txt` still has all results. Only organisms that share k-mers with the sample are reported, so this needs a very large database and `--show_all` or `--keep_raw`.

#### Parameters

The `--significance` parameter is basically akin to your confidence level: how sure do you want to be that the organism is present? Higher leads to more false negatives, lower leads to more false positives. 

The `--min_coverage_list` parameter dictates a list of `min_coverage` which indicates what percentage (value in `[0,1]`) of the distinct k-mers (think: whole genome) must have been sequenced and present in my sample to qualify as that organism as being "present." Setting this to 1 is usually safe, but if you have a very low coverage sample, you may want to lower this value. Setting it higher will lead to more false negatives, setting it lower will lead to more false positives (pretty rapidly).

| Parameter                               | Explanation                                                  |
| --------------------------------------- | ------------------------------------------------------------ |
| --json      | the path to a json file generated by `yacht train` (see above) |
| --significance                          | minimum probability of individual true negative, between 0 and 1 (default: 0.99) |
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

When we get the results folder from `yacht run`, you can run `yacht convert` to convert the YACHT result to other popular output formats (currently `cami`, `biom` and `graphplan` are supported).

__Note__: Before you run `yacht convert`, you need to prepare a TSV file `genome_to_taxid.tsv` containing two columns: genome ID (genome_id) and its corresponding taxid (taxid). An example can be found [here](demo/toy_genome_to_taxid.tsv). You need to prepare it according to the reference database genomes you used. The [pretrained databases](#pretrained-reference-databases) include it as `<prefix>_genome_to_taxid.tsv`. 

Then you are ready to run `yacht convert` with something like:
```bash
yacht convert --yacht_output_dir './results' --sheet_name 'min_coverage0.01' --genome_to_taxid 'genome_to_taxid.tsv' --mode 'cami' --sample_name 'MySample' --outfile_prefix 'cami_result' --outdir ./
```

#### Parameters

| Parameter         | Explanation                                                  |
| ----------------- | ------------------------------------------------------------ |
| --yacht_output_dir    | the path to the `results` folder generated by `yacht run` (containing `result.xlsx`) |
| --sheet_name      | specify which spreadsheet result you want to convert from     |
| --genome_to_taxid | the path to the location of `genome_to_taxid.tsv` you prepared |
| --mode            | specify to which output format you want to convert (e.g., 'cami', 'biom', 'graphplan')
| --sample_name     | A random name you would like to show in header of the cami file. Default: Sample1.' |
| --outfile_prefix  | the prefix of the output file. Default: result | 
| --outdir          | the path to output directory where the results will be generated |



