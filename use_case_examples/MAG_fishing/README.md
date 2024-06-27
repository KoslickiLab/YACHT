# YACHT for MAG Fishing

## Description of use-case-example

Metagenomic Assembled-Genomes (MAG) fishing is the process of reporting the assembled genomes within a metagenomic sample.

Metagenomics has been an important field in exploring the microbial communities of specific environments, especially for environments that contain unculturable microbes. However, there is a persistent underrepresentation of genomes challenging the production of a high-resolution of taxonomic profile. Consequently, many microbial communities are still understudied. Efforts have been made to increase the knowledge of these environments, such as the study we highlight here. One of the goals in the study by Banchi and colleagues (link to paper) was to unveil a more resolved taxonomic composition of marine sediments in the Venice Lagoon for further functional analyses of these microbial communities. The dataset from this study (NCBI accession: PRJNA924243) has 58 MAGS and serves as a use case example of using YACHT to resolve taxonomic composition.

According to their study, we should expect YACHT to report species from the phylum Proteobacteria, under the classes Alphaproteobacteria, Gammaproteobacteria, and Deltaprotwobacteria.

Banchi, E., Corre, E., Del Negro, P., Celussi, M., & Malfatti, F. (2024). Genome-resolved metagenomics of Venice Lagoon surface sediment bacteria reveals high biosynthetic potential and metabolic plasticity as successful strategies in an impacted environment. Marine Life Science & Technology, 6(1), 126-142.

## Install the following programs

**datasets**

More information please go to: [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

**YACHT**

More information please visit: [YACHT](https://github.com/KoslickiLab/YACHT).


## Download MAG samples

The following command was used to download the MAG sample of interest:

```
datasets download genome accession PRJNA924243
```

Downloading this MAG project will produce a directory of multiple pathways for each fasta file in this project but yacht wants one fasta file with running `yacht sketch sample` and yacht will let you know of this with the following error message:

```
ValueError: Please provide either one file for single-end reads or two files for paired-end reads.
```

A work around is running the following commands.

```
cd MAG_data
cp data/ncbi_dataset/data/GCA_02928*/*fna MAG_data/.
```

## yacht sketch MAG of interest

```
yacht sketch sample --infile MAG_sample.fna --kmer 31 --scaled 1000 --outfile sample.sig.zip
```

Be aware that `yacht sketch sample` will create a sketch sample with more than one signature, but `yacht run` wants a sample with one signature, so it will direct you to create merge signatures using `sourmash merge`. Please execute the following command:

```
sourmash sig merge sample.sig.zip -k 31 -o sample_merge.sig.zip
```

## Using yacht download pretrained_ref_db

Trying to download the pretrained data did not download anything

```
yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 31 --ani_thresh 0.9995 --outfolder ./
```

Runnning the following command, gave me everything?

```
yacht download default_ref_db --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder ./
```

It seems that the yacht downloand pretrained_ref_db doesn't show? It never really completed but completed once I ran yacht download default_ref_db

## yacht run
```
yacht run --json 'gtdb-rs214-reps.k31_0.9995_pretrained/gtdb-rs214-reps.k31_0.9995_config.json' --sample_file 'sample_merge.sig.zip' --num_threads 32 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out ./result.xlsx
```

