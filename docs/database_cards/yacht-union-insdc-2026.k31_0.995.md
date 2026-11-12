# YACHT pretrained database `yacht-union-insdc-2026.k31_0.995`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

Bacteria, archaea, fungi, protozoa and viruses, with INSDC/NCBI assembly accessions only: GTDB R232 species representatives, GenBank fungi and protozoa (April 2025, without the assemblies NCBI flags as contaminated), NCBI RefSeq viral genomes (January 2025) and the ICTV exemplar viruses (VMR MSL40). Every genome can be retrieved from NCBI by its accession; this is the easiest database to go from a result to the genome.

| | |
|---|---|
| genomes in the database | 230291 |
| genome sketches it was trained on | 256362 (26071 left out as equivalent to another genome at ANI >= 0.995) |
| ANI threshold | 0.995 |
| k-mer size / scaled | 31 / 1000 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.995`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `yacht-union-insdc-2026.k31_0.995_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | CC BY-SA 4.0 (as GTDB) |

## Use

```bash
yacht download pretrained_ref_db --name yacht-union-insdc-2026.k31 --ani_thresh 0.995 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht run --json yacht-union-insdc-2026.k31_0.995_pretrained/yacht-union-insdc-2026.k31_0.995_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 1000 also works (it is downsampled to 1000).

## Sources

- **GTDB R11-RS232 (release 232.0, 15 April 2026), species representative genomes (199,923)**
  - https://data.gtdb.ecogenomic.org/releases/release232/232.0/
  - files: genomic_files_reps/gtdb_genomes_reps_r232.tar.gz (md5 a63a586186871a429bf12ea5bc7be0db), bac120/ar53 metadata and taxonomy
  - license: CC BY-SA 4.0
  - please cite: Parks et al., GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Res. 50:D785 (2022), https://doi.org/10.1093/nar/gkab776
  - preparation: Sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000` (sourmash_plugin_branchwater), named `<accession> <ncbi_organism_name>` from the GTDB metadata (the naming of sourmash's GTDB databases); k=31 used for training.
- **GenBank fungal genomes, 8 April 2025 (21,135 genomes), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/ (genbank-20250408-fungi-k31.zip, lineages.fungi.csv)
  - license: NCBI GenBank (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI GenBank, https://www.ncbi.nlm.nih.gov/genbank/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: Used as published (k=31, scaled=1000), without the assemblies NCBI flags as contaminated (`excluded_from_refseq` contains "contaminated" in the NCBI GenBank assembly summary, current or historical, downloaded 30 October 2026), left out with `sourmash sig extract --picklist ...:exclude`; they are listed in `fungi.ncbi-contaminated.csv` in this zip.
- **GenBank protozoan genomes, 8 April 2025 (2,392 genomes), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/ (genbank-20250408-protozoa-k31.zip, lineages.protozoa.csv)
  - license: NCBI GenBank (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI GenBank, https://www.ncbi.nlm.nih.gov/genbank/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: Used as published (k=31, scaled=1000), without the assemblies NCBI flags as contaminated (`excluded_from_refseq` contains "contaminated" in the NCBI GenBank assembly summary, current or historical, downloaded 30 October 2026), left out with `sourmash sig extract --picklist ...:exclude`; they are listed in `protozoa.ncbi-contaminated.csv` in this zip.
- **NCBI RefSeq viral genomes of the NCBI Virus assemblies of January 2025 (16,280 of 230,973 assemblies, the `GCF_` accessions), prepared by the sourmash team**
  - https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/ncbi-viruses-2025.01/ (ncbi-viruses.dna.k=31.scaled=50.sig.zip, ncbi-viruses.lineages.csv)
  - license: NCBI (no NCBI restrictions on use); sourmash sketches: CC0
  - please cite: NCBI Virus, https://www.ncbi.nlm.nih.gov/labs/virus/; sourmash: Irber et al., J. Open Source Softw. 9:6830 (2024), https://doi.org/10.21105/joss.06830
  - preparation: The RefSeq (`GCF_`) sketches, selected with `sourmash sig extract --picklist` (k=31, scaled=50); for the unions downsampled to scaled 1000 with `sourmash sig downsample`.
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

GTDB R232 taxonomy for bacteria and archaea (`d__...`), NCBI taxonomy for eukaryotes and RefSeq viruses, ICTV MSL40 taxonomy for the VMR exemplars (in the same columns, `superkingdom` = `Viruses`). NCBI taxids for all genomes that have one (VMR exemplars: from the NCBI viral lineages of the same accession).

- `yacht-union-insdc-2026.k31_0.995_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `yacht-union-insdc-2026.k31_0.995_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 2798 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.995). For each taxonomic rank:
the taxa of the 256362 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 25780 left-out genomes that have a genome of this database at ANI >= 0.995:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 4 | 4 (100.0%) | 0.6% |
| phylum | 246 | 245 (99.6%) | 0.7% |
| class | 791 | 789 (99.7%) | 3.9% |
| order | 2,705 | 2,701 (99.9%) | 4.2% |
| family | 7,761 | 7,745 (99.8%) | 5.7% |
| genus | 43,110 | 42,953 (99.6%) | 7.9% |
| species | 228,965 | 224,831 (98.2%) | 31.3% |
| strain | 10,924 | 5,393 (49.4%) | 31.9% |

## Files

| file | content |
|---|---|
| `yacht-union-insdc-2026.k31_0.995_config.json`, `yacht-union-insdc-2026.k31_0.995_processed_manifest.tsv`, `yacht-union-insdc-2026.k31_0.995_index.yidx` | the database `yacht run` uses |
| `yacht-union-insdc-2026.k31_0.995_lineages.csv`, `yacht-union-insdc-2026.k31_0.995_genome_to_taxid.tsv` | taxonomy (see above) |
| `yacht-union-insdc-2026.k31_0.995_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.995; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.995 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.

## Notes

- **Viruses** are the NCBI RefSeq viral genomes and the ICTV exemplars (VMR MSL40): curated genomes only. The GenBank
  viral assemblies are not included, because some of them are collections of many viral genomes or contain their host's
  genome, and `yacht train` keeps the larger of two equivalent genomes. Viruses are included at scaled 1000, where a
  30 kb genome has about 30 hashes and many small viruses have none (those are left out: they can never be detected).
  For viruses, the scaled-50 viral databases (`ncbi-viruses-refseq-2025.01.k31-sc50`, `ictv-vmr-msl40.k31-sc50`, with
  samples sketched at `--scaled 50`) are more sensitive.
- A few ICTV exemplars are sketches of the host assembly that carries the virus (prophages, gene transfer agents, e.g.
  `GCA_000091085.2` "Chlamydia phage phiCPAR39", the 1.25 Mbp *Chlamydia pneumoniae* AR39 genome). Since `yacht train`
  keeps the larger of two equivalent genomes, such an exemplar can represent its host's GTDB genome (a handful of
  genomes; see `yacht-union-insdc-2026.k31_0.995_removed_genomes.tsv`).
- **Eukaryotes:** the GenBank fungal and protozoan assemblies NCBI flags as contaminated are left out (listed in
  `fungi.ncbi-contaminated.csv` and `protozoa.ncbi-contaminated.csv`). Other eukaryotic assemblies can still contain
  bacterial sequence, and since `yacht train` keeps the larger of two equivalent genomes, such an assembly can represent
  a bacterial genome: `yacht-union-insdc-2026.k31_0.995_removed_genomes.tsv` lists, for each left-out genome, the genome that represents it.
- **Prophages:** `yacht train` treats two genomes as equivalent when the smaller one is contained in the larger one at
  ANI >= 0.995. A phage whose genome is (as a prophage) part of a larger genome in the database is therefore represented by
  that genome (see `yacht-union-insdc-2026.k31_0.995_removed_genomes.tsv`).
- Taxonomies differ between the parts of the database: GTDB for bacteria and archaea, NCBI for eukaryotes and RefSeq
  viruses, ICTV MSL40 for the VMR exemplars.
