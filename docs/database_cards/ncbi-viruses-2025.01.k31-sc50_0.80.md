# YACHT pretrained database `ncbi-viruses-2025.01.k31-sc50_0.80`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

NCBI viral genome assemblies (January 2025), all GenBank and RefSeq assemblies, at scaled 50. See the note on GenBank assemblies that are not single viral genomes; `ncbi-viruses-refseq-2025.01.k31-sc50` holds the curated RefSeq genomes only.

| | |
|---|---|
| genomes in the database | 19218 |
| genome sketches it was trained on | 230967 (211749 left out as equivalent to another genome at ANI >= 0.80) |
| ANI threshold | 0.80 |
| k-mer size / scaled | 31 / 50 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.80`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `ncbi-viruses-2025.01.k31-sc50_0.80_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | NCBI data; sketches CC0 |

> **ANI 0.80 is a very coarse setting, much coarser than genus level.** At k = 31, ANI 0.80 corresponds to a
> k-mer containment of only 0.00099: two genomes are equivalent when about 1 in 1,000 hashes of the smaller one
> is shared (about 4 of the 4,000 hashes of a 4 Mbp genome at scaled 1000, a single hash of a 30 kb virus at scaled
> 50). Overlaps this small also occur between genomes of different families, orders and phyla. The table "How coarse
> is this database?" below shows how many taxa keep no genome of their own and how often a left-out genome is
> represented by a genome of another taxon. Use this database only for a coarse overview; for genus- or species-level
> results, use ANI 0.95 or higher.

## Use

```bash
yacht download pretrained_ref_db --name ncbi-viruses-2025.01.k31-sc50 --ani_thresh 0.80 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 50 --outfile sample.sig.zip
yacht run --json ncbi-viruses-2025.01.k31-sc50_0.80_pretrained/ncbi-viruses-2025.01.k31-sc50_0.80_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 50 also works (it is downsampled to 50).

## Sources

- **NCBI Virus genome assemblies (NCBI:txid10239), January 2025 (230,973 assemblies), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/ncbi-viruses-2025.01/ (ncbi-viruses.dna.k=31.scaled=50.sig.zip, ncbi-viruses.lineages.csv)
  - license: NCBI (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI Virus, https://www.ncbi.nlm.nih.gov/labs/virus/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: Used as published (k=31, scaled=50).

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCA_...` for GenBank, `GCF_...` for RefSeq), e.g. `GCF_000005845.2`.
How to retrieve a genome from its identifier: `datasets download genome accession <accession> --include genome` (NCBI Datasets), or https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

NCBI taxonomy and taxids (`ncbi-viruses.lineages.csv` of the sourmash team).

- `ncbi-viruses-2025.01.k31-sc50_0.80_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `ncbi-viruses-2025.01.k31-sc50_0.80_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 0 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.80). For each taxonomic rank:
the taxa of the 230967 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 96471 left-out genomes that have a genome of this database at ANI >= 0.80:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 1 | 1 (100.0%) | 0.0% |
| phylum | 17 | 17 (100.0%) | 0.1% |
| class | 41 | 41 (100.0%) | 0.1% |
| order | 73 | 72 (98.6%) | 0.1% |
| family | 259 | 248 (95.8%) | 0.2% |
| genus | 2,805 | 1,904 (67.9%) | 4.3% |
| species | 33,821 | 12,062 (35.7%) | 45.4% |
| strain | 61,169 | 4,213 (6.9%) | 78.2% |

## Files

| file | content |
|---|---|
| `ncbi-viruses-2025.01.k31-sc50_0.80_config.json`, `ncbi-viruses-2025.01.k31-sc50_0.80_processed_manifest.tsv`, `ncbi-viruses-2025.01.k31-sc50_0.80_index.yidx` | the database `yacht run` uses |
| `ncbi-viruses-2025.01.k31-sc50_0.80_lineages.csv`, `ncbi-viruses-2025.01.k31-sc50_0.80_genome_to_taxid.tsv` | taxonomy (see above) |
| `ncbi-viruses-2025.01.k31-sc50_0.80_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.80; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.80 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.

## Notes

- Viral genomes are small: at scaled 50, a 30 kb genome has about 600 hashes. The fewer hashes a genome
  has, the less evidence YACHT has for or against its presence.

- **Not every GenBank "viral" assembly is a single viral genome.** 92 assemblies of this set are larger than 1 Mbp,
  while the largest known viral genomes are about 2.5 Mbp: e.g. `GCA_024181985.1` "uncultured virus" (334 Mbp, a
  collection of metagenomic viral contigs) and `GCA_038430685.1` "Enterococcus phage Efs 1-1-1" (3 Mbp, with its
  host's genome). `yacht train` keeps the larger of two equivalent genomes, so such an assembly can represent many
  viruses (at ANI 0.95, `GCA_024181985.1` alone represents 5,042), and an assembly that contains a bacterial genome is
  reported present in samples with that bacterium. `ncbi-viruses-refseq-2025.01.k31-sc50` (the curated RefSeq genomes
  only) and `ictv-vmr-msl40.k31-sc50` avoid this.
