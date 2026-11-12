# YACHT pretrained database `gtdb-rs232-reps.k31_0.995`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

GTDB R232 species representative genomes (bacteria and archaea).

| | |
|---|---|
| genomes in the database | 199923 |
| genome sketches it was trained on | 199923 (0 left out as equivalent to another genome at ANI >= 0.995) |
| ANI threshold | 0.995 |
| k-mer size / scaled | 31 / 1000 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.995`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `gtdb-rs232-reps.k31_0.995_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | CC BY-SA 4.0 (as GTDB) |

## Use

```bash
yacht download pretrained_ref_db --name gtdb-rs232-reps.k31 --ani_thresh 0.995 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht run --json gtdb-rs232-reps.k31_0.995_pretrained/gtdb-rs232-reps.k31_0.995_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 1000 also works (it is downsampled to 1000).

## Sources

- **GTDB R11-RS232 (release 232.0, 15 April 2026), species representative genomes (199,923)**
  - https://data.gtdb.ecogenomic.org/releases/release232/232.0/
  - files: genomic_files_reps/gtdb_genomes_reps_r232.tar.gz (md5 a63a586186871a429bf12ea5bc7be0db), bac120/ar53 metadata and taxonomy
  - license: CC BY-SA 4.0
  - please cite: Parks et al., GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Res. 50:D785 (2022), https://doi.org/10.1093/nar/gkab776
  - preparation: Sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000` (sourmash_plugin_branchwater), named `<accession> <ncbi_organism_name>` from the GTDB metadata (the naming of sourmash's GTDB databases); k=31 used for training.

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCA_...` for GenBank, `GCF_...` for RefSeq), e.g. `GCF_000005845.2`.
How to retrieve a genome from its identifier: `datasets download genome accession <accession> --include genome` (NCBI Datasets), or https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/. Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

GTDB R232 taxonomy (domain to species, `d__...;s__...`), from the GTDB metadata; NCBI taxids from the GTDB metadata (`ncbi_taxid`).

- `gtdb-rs232-reps.k31_0.995_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `gtdb-rs232-reps.k31_0.995_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 0 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.995). For each taxonomic rank:
the taxa of the 199923 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 0 left-out genomes that have a genome of this database at ANI >= 0.995:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 2 | 2 (100.0%) | - |
| phylum | 197 | 197 (100.0%) | - |
| class | 645 | 645 (100.0%) | - |
| order | 2,351 | 2,351 (100.0%) | - |
| family | 6,747 | 6,747 (100.0%) | - |
| genus | 37,503 | 37,503 (100.0%) | - |
| species | 199,923 | 199,923 (100.0%) | - |

## Files

| file | content |
|---|---|
| `gtdb-rs232-reps.k31_0.995_config.json`, `gtdb-rs232-reps.k31_0.995_processed_manifest.tsv`, `gtdb-rs232-reps.k31_0.995_index.yidx` | the database `yacht run` uses |
| `gtdb-rs232-reps.k31_0.995_lineages.csv`, `gtdb-rs232-reps.k31_0.995_genome_to_taxid.tsv` | taxonomy (see above) |
| `gtdb-rs232-reps.k31_0.995_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.995; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.995 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.
