# YACHT pretrained database `genbank-2025.04-fungi-k31_0.9995`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

GenBank fungal genomes (April 2025).

| | |
|---|---|
| genomes in the database | 15564 |
| genome sketches it was trained on | 21070 (5506 left out as equivalent to another genome at ANI >= 0.9995) |
| ANI threshold | 0.9995 |
| k-mer size / scaled | 31 / 1000 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.9995`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `genbank-2025.04-fungi-k31_0.9995_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | NCBI data; sketches CC0 |

## Use

```bash
yacht download pretrained_ref_db --name genbank-2025.04-fungi-k31 --ani_thresh 0.9995 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht run --json genbank-2025.04-fungi-k31_0.9995_pretrained/genbank-2025.04-fungi-k31_0.9995_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 1000 also works (it is downsampled to 1000).

## Sources

- **GenBank fungal genomes, 8 April 2025 (21,135 genomes), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/ (genbank-20250408-fungi-k31.zip, lineages.fungi.csv)
  - license: NCBI GenBank (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI GenBank, https://www.ncbi.nlm.nih.gov/genbank/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: Used as published (k=31, scaled=1000), without the assemblies NCBI flags as contaminated (`excluded_from_refseq` contains "contaminated" in the NCBI GenBank assembly summary, current or historical, downloaded 30 October 2026), left out with `sourmash sig extract --picklist ...:exclude`; they are listed in `fungi.ncbi-contaminated.csv` in this zip.

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCA_...`), e.g. `GCF_000005845.2`.
How to retrieve a genome from its identifier: `datasets download genome accession <accession> --include genome` (NCBI Datasets), or https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

NCBI taxonomy and taxids (`lineages.fungi.csv` of the sourmash team).

- `genbank-2025.04-fungi-k31_0.9995_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `genbank-2025.04-fungi-k31_0.9995_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 0 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.9995). For each taxonomic rank:
the taxa of the 21070 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 5441 left-out genomes that have a genome of this database at ANI >= 0.9995:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 1 | 1 (100.0%) | 0.0% |
| phylum | 10 | 10 (100.0%) | 0.0% |
| class | 55 | 55 (100.0%) | 0.2% |
| order | 166 | 166 (100.0%) | 0.6% |
| family | 490 | 490 (100.0%) | 0.8% |
| genus | 1,448 | 1,442 (99.6%) | 1.1% |
| species | 6,405 | 6,218 (97.1%) | 8.3% |
| strain | 1,041 | 818 (78.6%) | 8.6% |

## Files

| file | content |
|---|---|
| `genbank-2025.04-fungi-k31_0.9995_config.json`, `genbank-2025.04-fungi-k31_0.9995_processed_manifest.tsv`, `genbank-2025.04-fungi-k31_0.9995_index.yidx` | the database `yacht run` uses |
| `genbank-2025.04-fungi-k31_0.9995_lineages.csv`, `genbank-2025.04-fungi-k31_0.9995_genome_to_taxid.tsv` | taxonomy (see above) |
| `genbank-2025.04-fungi-k31_0.9995_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.9995; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.9995 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.

## Notes

- The assemblies NCBI flags as contaminated are left out (see Sources). Other GenBank eukaryotic assemblies can still
  contain contaminating (e.g. bacterial) sequence. With this eukaryote-only
  database, a genome reported present with many matches in a bacteria-rich sample may be detected through such
  sequence. In the union databases, k-mers shared with the bacterial genomes are not exclusive to the eukaryotic genome,
  so they do not count as evidence for it.
