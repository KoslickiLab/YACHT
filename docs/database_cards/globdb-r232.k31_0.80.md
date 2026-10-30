# YACHT pretrained database `globdb-r232.k31_0.80`

**For YACHT 2.0.0 or later.** YACHT 1.x cannot use this database.

GlobDB r232: the most comprehensive species-level set of bacterial and archaeal genomes (GTDB R232 plus MAG catalogues).

| | |
|---|---|
| genomes in the database | 17558 |
| genome sketches it was trained on | 346233 (328675 left out as equivalent to another genome at ANI >= 0.80) |
| ANI threshold | 0.80 |
| k-mer size / scaled | 31 / 1000 |
| trained with | YACHT 2.0.0 (`yacht train --ksize 31 --ani_thresh 0.80`) |
| requires | YACHT 2.0.0 or later (`min_yacht_version` in `globdb-r232.k31_0.80_config.json`; YACHT 1.x needs the genome sketches, which are not included) |
| license | CC BY-SA 4.0 (as GlobDB and GTDB) |

> **ANI 0.80 is a very coarse setting, much coarser than genus level.** At k = 31, ANI 0.80 corresponds to a
> k-mer containment of only 0.00099: two genomes are equivalent when about 1 in 1,000 hashes of the smaller one
> is shared (about 4 of the 4,000 hashes of a 4 Mbp genome at scaled 1000, a single hash of a 30 kb virus at scaled
> 50). Overlaps this small also occur between genomes of different families, orders and phyla. The table "How coarse
> is this database?" below shows how many taxa keep no genome of their own and how often a left-out genome is
> represented by a genome of another taxon. Use this database only for a coarse overview; for genus- or species-level
> results, use ANI 0.95 or higher.

## Use

```bash
yacht download pretrained_ref_db --name globdb-r232.k31 --ani_thresh 0.80 --outfolder ./
yacht sketch sample --infile reads_1.fq.gz reads_2.fq.gz --kmer 31 --scaled 1000 --outfile sample.sig.zip
yacht run --json globdb-r232.k31_0.80_pretrained/globdb-r232.k31_0.80_config.json --sample_file sample.sig.zip --outdir ./
```

A sample sketched with a smaller `--scaled` than 1000 also works (it is downsampled to 1000).

## Sources

- **GTDB R11-RS232 (release 232.0, 15 April 2026), species representative genomes (199,923)**
  - https://data.gtdb.ecogenomic.org/releases/release232/232.0/
  - files: genomic_files_reps/gtdb_genomes_reps_r232.tar.gz (md5 a63a586186871a429bf12ea5bc7be0db), bac120/ar53 metadata and taxonomy
  - license: CC BY-SA 4.0
  - please cite: Parks et al., GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy, Nucleic Acids Res. 50:D785 (2022), https://doi.org/10.1093/nar/gkab776
  - preparation: Sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000` (sourmash_plugin_branchwater), named `<accession> <ncbi_organism_name>` from the GTDB metadata (the naming of sourmash's GTDB databases); k=31 used for training.
- **GlobDB release 232 (June 2026): 346,233 species-representative genomes, the GTDB R232 species representatives plus species-level MAGs from 25 further catalogues**
  - https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/
  - files: globdb_r232_genome_fasta.tar.gz (md5 c95ef3b3278b9f035e2089bd729b0b48), globdb_r232_taxonomy.tsv.gz, globdb_r232_dictionaries.tar.gz
  - license: CC BY-SA 4.0 (GlobDB propagates the licenses of the underlying data sources)
  - please cite: Speth et al., GlobDB: a comprehensive species-dereplicated microbial genome resource, Bioinformatics Advances (2025), https://doi.org/10.1093/bioadv/vbaf280, and the source catalogues listed at https://globdb.org/citations
  - preparation: The GTDB R232 part reuses the GTDB R232 sketches above (same genomes, versioned accessions from the GTDB metadata); the other genomes were sketched with `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000`, named `<GlobDB ID> <GlobDB species>`.

## Genome identifiers

The `organism_name` of each genome in YACHT's output is its sketch name. Its first word is the genome's identifier:
a versioned INSDC/NCBI assembly accession (`GCA_...`/`GCF_...`) for the GTDB R232 genomes, and the GlobDB genome ID (e.g. `MOTU40_133241`, `SPIREOTU_03021551`, `MGYG000535630`; the prefix names the source catalogue) for the other genomes. GlobDB IDs are stable within a GlobDB release; `globdb_r232_dictionaries.tar.gz` maps them to the IDs of the source catalogues.
How to retrieve a genome from its identifier: for an accession, NCBI Datasets (`datasets download genome accession <accession>`); for a GlobDB ID, the genome is `globdb_r232_genome_fasta/<chunk>/<GlobDB ID>.fa.gz` in `globdb_r232_genome_fasta.tar.gz` of the GlobDB r232 release (it can be extracted from the download stream without storing the whole archive). Details: https://github.com/KoslickiLab/YACHT/blob/main/docs/retrieving_genomes.md

## Taxonomy

GlobDB r232 taxonomy (GTDB R232 taxonomy, assigned with GTDB-Tk for the non-GTDB genomes). NCBI taxids only for the GTDB R232 genomes (from the GTDB metadata).

- `globdb-r232.k31_0.80_lineages.csv`: lineage of every genome in the database (sourmash lineage format; `ident` = identifier).
  0 genomes have no lineage.
- `globdb-r232.k31_0.80_genome_to_taxid.tsv`: NCBI taxid of every genome that has one (`genome_id` = `organism_name`), the input
  `yacht convert --genome_to_taxid` needs for CAMI/BIOM output. 5997 genomes have no NCBI taxid.

## How coarse is this database?

`yacht train` keeps one genome of each group of genomes it finds equivalent (ANI >= 0.80). For each taxonomic rank:
the taxa of the 346233 genomes it was trained on, how many of them still have a genome in this database, and how
often the closest database genome of a left-out genome belongs to another taxon at that rank (or above), counted over
the 180433 left-out genomes that have a genome of this database at ANI >= 0.80:

| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |
|---|---:|---:|---:|
| superkingdom | 2 | 2 (100.0%) | 0.9% |
| phylum | 223 | 204 (91.5%) | 15.4% |
| class | 760 | 629 (82.8%) | 18.7% |
| order | 3,058 | 1,904 (62.3%) | 29.5% |
| family | 12,090 | 5,446 (45.0%) | 38.8% |
| genus | 75,469 | 14,905 (19.7%) | 71.1% |
| species | 346,233 | 17,558 (5.1%) | 100.0% |

## Files

| file | content |
|---|---|
| `globdb-r232.k31_0.80_config.json`, `globdb-r232.k31_0.80_processed_manifest.tsv`, `globdb-r232.k31_0.80_index.yidx` | the database `yacht run` uses |
| `globdb-r232.k31_0.80_lineages.csv`, `globdb-r232.k31_0.80_genome_to_taxid.tsv` | taxonomy (see above) |
| `globdb-r232.k31_0.80_removed_genomes.tsv` | each genome `yacht train` left out, the closest genome of this database it is equivalent to (ANI >= 0.80; highest containment, with the ANI estimate containment^(1/k)), and how many genomes of the database it is equivalent to. `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is equivalent to a larger genome that has not been left out; that larger genome can be left out later in turn. A genome left out through such a chain can have no genome of this database at ANI >= 0.80 (empty closest genome, 0 equivalent genomes). |

The reference genome sketches themselves are not included (`yacht run` only needs the index); see Sources.
