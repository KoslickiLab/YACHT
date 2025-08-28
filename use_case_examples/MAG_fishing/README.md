
# MAG Fishing

## Introduction

YACHT can be used to search or "fish" for a metagenome-assembled genome (MAG) of interest in a sample. We take a closer look of the results produced for the use case example "Low Abundance Samples", and find that *Bacteroidales bacterium* was present in sample SRR32008482 using a *k*-size of 51 and ANI of 0.95. In Table 1, we observe that the accepted threshold was to match at least 75 exclusive *k*-mer matches to the reference. Sample SRR32008482 exceeded this threshold, matching 158 exclusive *k*-mers.. Please refer to `result_k51_ani0.95_SRR32008482.xlsx` for more YACHT results of this sample. In the following, we are curious to find out whether another metagenomic sample, SRR32008483, has *B. bacterium* present in its microbial community.

### Table 1: B. bacterium detected in Sample SRR32008482 using k-size of 51 and ani=0.95
| Sample      | Minimum Coverage | Organism Name                                         | Exclusive *k*-mers to Genome | Coverage of Exclusive *k*-mers | num_matches | Acceptance Threshold |
|-------------|------------------|-------------------------------------------------------|----------------------------|------------------------------|-------------|----------------------|
| SRR32008482 | 0.5              | GCA_017506175.1 Bacteroidales bacterium, ASM1750617v1 | 2464                       | 1232                         | 158         | 75                   |


Make sure you have the following dependencies to run this use case example:

- YACHT
- fastq-dump
- datasets
- matplotlib_venn

## Obtain datasets

Throughout this use case example, we will use two sample datasets, SRR32008482 and SRR32008483, to evaluate how results may differ in samples when MAG fishing. Note that these two samples have different sequence coverages where one has a sequence coverage of 44Mb (SRR32008482) and the other has a sequence coverage of 4.4Mb (SRR32008483). Additionally, we will use only one genome of interest as our MAG reference dataset to train with YACHT, *B. bacterium* (Accession: GCA_017506175).

Please run the following script to download these samples and the MAG reference dataset.

    bash data_download.sh

## Can we find the MAG, *Bacteroidales bacterium*, in both of these samples?

### Sketch samples of interest

We start off our analysis by sketching samples using `yacht sketch sample`. 

#### Sample SRR32008482:
    yacht sketch sample --infile data/SRR32008482.fastq --kmer 51 --scaled 10 --outfile SRR32008482.k51.sample.sig.zip
#### Sample SRR32008483:
    yacht sketch sample --infile data/SRR32008483.fastq --kmer 51 --scaled 10 --outfile SRR32008483.k51.sample.sig.zip

<!--

## For internal purposes, I want to identify whether these two share a species

## Sketch samples of interest

We start off our analysis by sketching samples using `yacht sketch sample`. 

### Sample SRR32008482:
    yacht sketch sample --infile data/SRR32008482.fastq --kmer 51 --scaled 1000 --outfile SRR32008482.k51.sample.sig.zip
### Sample SRR32008483:
    yacht sketch sample --infile data/SRR32008483.fastq --kmer 51 --scaled 1000 --outfile SRR32008483.k51.sample.sig.zip


I downloaded pre-trained reference signature using `yacht download`.

    yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 51 --ani_thresh 0.95 --outfolder ./

Now I run `yacht run` for both samples

    nohup yacht run --json 'gtdb-rs214-reps.k51_0.95_pretrained/gtdb-rs214-reps.k51_0.95_config.json' --sample_file 'SRR32008482.k51.sample.sig.zip' --num_threads 32 --keep_raw --significance 0.95 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result_k51_ani0.95_SRR32008482.xlsx > run.k51.SRR32008482.log 2>&1 &

    nohup yacht run --json 'gtdb-rs214-reps.k51_0.95_pretrained/gtdb-rs214-reps.k51_0.95_config.json' --sample_file 'SRR32008483.k51.sample.sig.zip' --num_threads 32 --keep_raw --significance 0.95 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result_k51_ani0.95_SRR32008483.xlsx > run.k51.SRR32008483.log 2>&1 &

Do these share any genomes?

Using a k-size of 51, we look at results for a minimum coverage of 0.05. Both of these samples share 21 species, where Sample SRR32008483 has 9 unique species and Sample SRR32008483 has 200 unique species.

    python venn.py

![Alt text](venn.png)

We have a lot of candidates. I chose GCA_017506175.1_ASM1750617v1_genomic.fna

-->

### Sketch the MAG of interest, *B. bacterium*

If you completed the "Low abundance sample" use case example, you may recall that we were able to download pre-trained reference datasets. However, in this use case example, we are only interested in retrieving or "fishing" a single reference MAG, rather than an entire reference database like GTDB. Therefore, we can use a smaller scale factor (i.e., 10) and still ensure computational efficiency preserving the *k*-mer set preventing any error that YACHT cannot match a *k*-mer from the MAG to our samples. 

Additionally, this reference will be a sketch using "--singleton" meaning that each entry within the FASTA/Q file will have a unique signature. Specifically, every entry in this reference file is a scaffold indicating that this reference has not yet been appropriately annotated.

Note that when sketching the sample or reference, that they share the same scale factor.

#### Similar to sketching a sample, we use `yacht sketch ref` to sketch our MAG of interest

    yacht sketch ref --infile data/ncbi_dataset/data/GCA_017506175.1/GCA_017506175.1_ASM1750617v1_genomic.fna --kmer 51 --scaled 10 --outfile MAG_ref_database.k51.sig.zip

### Using `yacht train`, we train our MAG signature
    yacht train --ref_file MAG_ref_database.k51.sig.zip --ksize 51 --num_threads 32 --ani_thresh 0.95 --prefix 'MAG_ref_database_k51_ani0.95' --outdir ./ 

### Identify whether *B. bacterium* is present or absent using `yacht run`

#### Initially, we were interested in whether *B. bacterium* was present in Sample SRR32008483.
    yacht run --json 'MAG_ref_database_k51_ani0.95_config.json' --sample_file 'SRR32008483.k51.sample.sig.zip' --num_threads 32 --keep_raw --significance 0.95 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result_SRR32008483.k51.xlsx

#### Because signatures were produced for each scaffold in the MAG reference database, we are also interested in evaluating the scaffolds present in the old sample of interest, Sample SRR32008482.
    yacht run --json 'MAG_ref_database_k51_ani0.95_config.json' --sample_file 'SRR32008482.k51.sample.sig.zip' --num_threads 32 --keep_raw --significance 0.95 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result_SRR32008482.k51.xlsx


## Do these samples share contigs of *B. bacterium*?

Run following script to produce and evaluate figure.

    python venn.py

Yes, we are able to identify contigs of *B. bacterium* within these samples. Both samples share 21 contigs of *B. bacterium*, where Sample SRR32008482 has 41 unique contigs at a minimal coverage of 0.5 and at a lower coverage Sample SRR32008483 has 2 unique contigs. Note, that the difference in scaffolds present can be an indication of difference in sequence coverage of the samples.

![Venn Diagram](venn.png)

## Reference

Degregori, Samuel, Melissa B. Manus, Evan B. Qu, Calen P. Mendall, Jacob S. Baker, Lydia M. Hopper, Katherine R. Amato, and Tami D. Lieberman. "The microbiome of the human facial skin is unique compared to that of other hominids." bioRxiv (2025): 2025-01.
