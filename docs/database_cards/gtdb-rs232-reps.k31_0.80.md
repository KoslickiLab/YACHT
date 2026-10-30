# YACHT pretrained database `gtdb-rs232-reps.k31_0.80`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

GTDB R232 species representative genomes (bacteria and archaea).

| | |
|---|---|
| genomes in the database | 14568 |
| genome sketches it was trained on | 199923 (185355 left out as equivalent to another genome at ANI >= 0.80) |
| ANI threshold | 0.80 |
| k-mer size / scaled | 31 / 1000 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.80`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `gtdb-rs232-reps.k31_0.80_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | CC BY-SA 4.0 (as GTDB) |

> **ANI 0.80 is a very coarse setting, much coarser than genus level.** At k = 31, ANI 0.80 corresponds to a
> k-mer containment of only 0.00099: two genomes are equivalent when about 1 in 1,000 hashes of the smaller one
> is shared (about 4 of the 4,000 hashes of a 4 Mbp genome at scaled 1000, a single hash of a 30 kb virus at scaled
> 50). Overlaps this small also occur between genomes of different families, orders and phyla. The table "How coarse
> is this database?" below shows how many taxa keep no genome of their own and how often a left-out genome is
> represented by a genome of another taxon. Use this database only for a coarse overview; for genus- or species-level
> results, use ANI 0.95 or higher.

## Use

```bash
yacht download pretrained_ref_db --name gtdb-rs232-reps.k31 --ani_thresh 0.80 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht run --json gtdb-rs232-reps.k31_0.80_pretrained/gtdb-rs232-reps.k31_0.80_config.json --sample_file sample.sig.zip --outdir ./
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

- `gtdb-rs232-reps.k31_0.80_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `gtdb-rs232-reps.k31_0.80_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 0 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.80). For each taxonomic rank:
the taxa of the 199923 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 104336 left-out genomes that have a genome of this database at ANI >= 0.80:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 2 | 2 (100.0%) | 0.8% |
| phylum | 197 | 182 (92.4%) | 13.2% |
| class | 645 | 556 (86.2%) | 15.8% |
| order | 2,351 | 1,625 (69.1%) | 27.3% |
| family | 6,747 | 3,976 (58.9%) | 36.8% |
| genus | 37,503 | 11,444 (30.5%) | 65.7% |
| species | 199,923 | 14,568 (7.3%) | 100.0% |

## Files

| file | content |
|---|---|
| `gtdb-rs232-reps.k31_0.80_config.json`, `gtdb-rs232-reps.k31_0.80_processed_manifest.tsv`, `gtdb-rs232-reps.k31_0.80_index.yidx` | the database `yacht run` uses |
| `gtdb-rs232-reps.k31_0.80_lineages.csv`, `gtdb-rs232-reps.k31_0.80_genome_to_taxid.tsv` | taxonomy (see above) |
| `gtdb-rs232-reps.k31_0.80_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.80; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.80 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.
