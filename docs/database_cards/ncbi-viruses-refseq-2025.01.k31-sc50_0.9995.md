# YACHT pretrained database `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

NCBI RefSeq viral genomes (January 2025), at scaled 50: the curated subset of the NCBI viral assemblies.

| | |
|---|---|
| genomes in the database | 16081 |
| genome sketches it was trained on | 16279 (198 left out as equivalent to another genome at ANI >= 0.9995) |
| ANI threshold | 0.9995 |
| k-mer size / scaled | 31 / 50 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.9995`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | NCBI data; sketches CC0 |

## Use

```bash
yacht download pretrained_ref_db --name ncbi-viruses-refseq-2025.01.k31-sc50 --ani_thresh 0.9995 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 50 --outfile sample.sig.zip
yacht run --json ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_pretrained/ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 50 also works (it is downsampled to 50).

## Sources

- **NCBI RefSeq viral genomes of the NCBI Virus assemblies of January 2025 (16,280 of 230,973 assemblies, the `GCF_` accessions), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/ncbi-viruses-2025.01/ (ncbi-viruses.dna.k=31.scaled=50.sig.zip, ncbi-viruses.lineages.csv)
  - license: NCBI (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI Virus, https://www.ncbi.nlm.nih.gov/labs/virus/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: The RefSeq (`GCF_`) sketches, selected with `sourmash sig extract --picklist` (k=31, scaled=50); for the unions downsampled to scaled 1000 with `sourmash sig downsample`.

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCF_...`, RefSeq), e.g. `GCF_000005845.2`.
How to retrieve a genome from its identifier: `datasets download genome accession <accession> --include genome` (NCBI Datasets), or https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

NCBI taxonomy and taxids (`ncbi-viruses.lineages.csv` of the sourmash team).

- `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 0 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.9995). For each taxonomic rank:
the taxa of the 16279 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 193 left-out genomes that have a genome of this database at ANI >= 0.9995:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 1 | 1 (100.0%) | 0.0% |
| phylum | 17 | 17 (100.0%) | 0.0% |
| class | 40 | 40 (100.0%) | 0.0% |
| order | 71 | 71 (100.0%) | 0.0% |
| family | 251 | 251 (100.0%) | 0.0% |
| genus | 2,759 | 2,758 (100.0%) | 0.0% |
| species | 14,129 | 14,074 (99.6%) | 31.6% |
| strain | 9,920 | 9,823 (99.0%) | 61.1% |

## Files

| file | content |
|---|---|
| `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_config.json`, `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_processed_manifest.tsv`, `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_index.yidx` | the database `yacht run` uses |
| `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_lineages.csv`, `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_genome_to_taxid.tsv` | taxonomy (see above) |
| `ncbi-viruses-refseq-2025.01.k31-sc50_0.9995_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.9995; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.9995 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.

## Notes

- Viral genomes are small: at scaled 50, a 30 kb genome has about 600 hashes. The fewer hashes a genome
  has, the less evidence YACHT has for or against its presence.
