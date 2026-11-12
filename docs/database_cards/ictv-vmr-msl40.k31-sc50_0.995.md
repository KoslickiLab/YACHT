# YACHT pretrained database `ictv-vmr-msl40.k31-sc50_0.995`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

ICTV exemplar virus genomes (Virus Metadata Resource MSL40), at scaled 50.

| | |
|---|---|
| genomes in the database | 17530 |
| genome sketches it was trained on | 17764 (234 left out as equivalent to another genome at ANI >= 0.995) |
| ANI threshold | 0.995 |
| k-mer size / scaled | 31 / 50 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.995`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `ictv-vmr-msl40.k31-sc50_0.995_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | CC BY 4.0 (as the ICTV VMR) |

## Use

```bash
yacht download pretrained_ref_db --name ictv-vmr-msl40.k31-sc50 --ani_thresh 0.995 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 50 --outfile sample.sig.zip
yacht run --json ictv-vmr-msl40.k31-sc50_0.995_pretrained/ictv-vmr-msl40.k31-sc50_0.995_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 50 also works (it is downsampled to 50).

## Sources

- **ICTV Virus Metadata Resource, MSL40 v1 (2025): 17,765 exemplar virus genomes, prepared by the sourmash team (ICTV has since released MSL41, July 2026)**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/2025-ictv-vmr-msl40/ (vmr_MSL40_v1.dna-k31.sig.zip, vmr_MSL40_v1.lineages.csv)
  - license: ICTV VMR: CC BY 4.0 (https://doi.org/10.5281/zenodo.15042310); sourmash sketches: CC0
  - please cite: ICTV Virus Metadata Resource, https://ictv.global/vmr, https://doi.org/10.5281/zenodo.15042310
  - preparation: Used as published (k=31, scaled=50); for the unions downsampled to scaled 1000 with `sourmash sig downsample`.

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCA_...` for GenBank, `GCF_...` for RefSeq), e.g. `GCF_000005845.2`.
How to retrieve a genome from its identifier: `datasets download genome accession <accession> --include genome` (NCBI Datasets), or https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

ICTV taxonomy (realm to species, MSL40; `vmr_MSL40_v1.lineages.csv` of the sourmash team). NCBI taxids of the assemblies from `ncbi-viruses.lineages.csv` (NCBI Virus, January 2025, same versioned accession).

- `ictv-vmr-msl40.k31-sc50_0.995_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `ictv-vmr-msl40.k31-sc50_0.995_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 2909 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.995). For each taxonomic rank:
the taxa of the 17764 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 231 left-out genomes that have a genome of this database at ANI >= 0.995:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| realm | 7 | 7 (100.0%) | 0.0% |
| kingdom | 11 | 11 (100.0%) | 0.0% |
| phylum | 22 | 22 (100.0%) | 0.0% |
| subphylum | 4 | 4 (100.0%) | 0.0% |
| class | 49 | 49 (100.0%) | 0.0% |
| order | 93 | 93 (100.0%) | 0.0% |
| suborder | 12 | 12 (100.0%) | 0.0% |
| family | 367 | 367 (100.0%) | 0.0% |
| subfamily | 213 | 213 (100.0%) | 0.0% |
| genus | 3,765 | 3,765 (100.0%) | 0.0% |
| subgenus | 86 | 86 (100.0%) | 0.0% |
| species | 16,089 | 15,982 (99.3%) | 46.8% |

## Files

| file | content |
|---|---|
| `ictv-vmr-msl40.k31-sc50_0.995_config.json`, `ictv-vmr-msl40.k31-sc50_0.995_processed_manifest.tsv`, `ictv-vmr-msl40.k31-sc50_0.995_index.yidx` | the database `yacht run` uses |
| `ictv-vmr-msl40.k31-sc50_0.995_lineages.csv`, `ictv-vmr-msl40.k31-sc50_0.995_genome_to_taxid.tsv` | taxonomy (see above) |
| `ictv-vmr-msl40.k31-sc50_0.995_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.995; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.995 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.

## Notes

- Viral genomes are small: at scaled 50, a 30 kb genome has about 600 hashes. The fewer hashes a genome
  has, the less evidence YACHT has for or against its presence.

- **A few ICTV exemplars are sketched with their host's genome.** For some exemplars (prophages, gene transfer
  agents) the VMR points to the GenBank assembly of the host that carries the virus, e.g. `GCA_000091085.2` "Chlamydia
  phage phiCPAR39" (the 1.25 Mbp *Chlamydia pneumoniae* AR39 genome) and `GCA_002951655.1` "Escherichia coli strain
  13P477T plasmid p1" (5.1 Mbp). Such an exemplar is reported present in samples that contain its host.