# YACHT 2.0 pretrained reference databases

All databases are trained with `yacht train` (YACHT 2.0) at k=31 for four ANI thresholds (0.80, 0.95, 0.995,
0.9995), and published in the [YACHT Zenodo community](https://zenodo.org/communities/yacht). List and download them
with:

```bash
yacht download pretrained_ref_db --list
yacht download pretrained_ref_db --name <name> --ani_thresh <ANI> --outfolder ./
```

`yacht download pretrained_ref_db` looks the databases up on Zenodo when it runs, so databases added later can be
downloaded without updating YACHT.

| name | content | genomes before training | identifiers | taxonomy | scaled |
|---|---|---:|---|---|---:|
| `gtdb-rs232-reps.k31` | GTDB R232 species representatives (bacteria, archaea) | 199,923 | versioned GCA/GCF accessions | GTDB R232 | 1000 |
| `globdb-r232.k31` | GlobDB r232: GTDB R232 species representatives + species-level MAGs of 25 catalogues | 346,233 | GCA/GCF accessions + GlobDB IDs | GlobDB r232 (GTDB R232) | 1000 |
| `genbank-2025.04-fungi-k31` | GenBank fungal genomes, April 2025 | 21,135 | GCA accessions | NCBI | 1000 |
| `genbank-2025.04-protozoa-k31` | GenBank protozoan genomes, April 2025 | 2,392 | GCA accessions | NCBI | 1000 |
| `ictv-vmr-msl40.k31-sc50` | ICTV exemplar viruses (VMR MSL40) | 17,765 | GCA accessions | ICTV MSL40 | 50 |
| `ncbi-viruses-2025.01.k31-sc50` | NCBI viral genome assemblies, January 2025 | 230,973 | GCA/GCF accessions | NCBI | 50 |
| `yacht-union-insdc-2026.k31` | GTDB R232 reps + GenBank fungi + protozoa + NCBI viruses | ~455,000 | GCA/GCF accessions only | GTDB (bacteria, archaea), NCBI (eukaryotes, viruses) | 1000 |
| `yacht-union-globdb-2026.k31` | GlobDB r232 + GenBank fungi + protozoa + NCBI viruses | ~600,000 | GCA/GCF accessions + GlobDB IDs | GlobDB (bacteria, archaea), NCBI (eukaryotes, viruses) | 1000 |

The numbers of genomes in each database after training are in the database's `README.md` (and in `--list`).

## Which database?

- **Bacteria and archaea, easiest to follow up:** `gtdb-rs232-reps.k31`. Every genome has an NCBI accession.
- **Bacteria and archaea, most comprehensive:** `globdb-r232.k31`. It adds about 146,000 species that are only known
  from metagenome-assembled genomes. Their genomes are retrieved from the GlobDB release instead of NCBI (see
  [retrieving genomes](retrieving_genomes.md)).
- **Everything (bacteria, archaea, fungi, protists, viruses):** `yacht-union-insdc-2026.k31` (every genome retrievable
  from NCBI by its accession; the easiest to follow up) or `yacht-union-globdb-2026.k31` (most comprehensive; GlobDB
  genomes are retrieved from the GlobDB release).
- **Viruses:** the scaled-50 databases `ictv-vmr-msl40.k31-sc50` (curated exemplars with ICTV taxonomy) and
  `ncbi-viruses-2025.01.k31-sc50` (all NCBI viral assemblies). Sketch the sample with `--scaled 50`; the same sketch
  also works with the scaled-1000 databases, since YACHT downsamples it.

## Identifiers and taxonomy (FAIR)

No identifiers were invented: every genome keeps the identifier of the resource it comes from.

- The first word of each genome's name (`organism_name` in YACHT's output) is its identifier: a versioned INSDC/NCBI
  assembly accession (`GCA_...`, `GCF_...`), or for the MAGs of GlobDB, the GlobDB r232 genome ID. Names are unique in
  every database.
- `<prefix>_lineages.csv` (sourmash lineage format) gives the taxonomy of every genome in the database:
  GTDB R232 for the GTDB genomes, the GlobDB r232 taxonomy for GlobDB, NCBI taxonomy for the GenBank and NCBI viral
  genomes, ICTV MSL40 for the VMR. In the unions, the column `superkingdom` holds `d__Bacteria`/`d__Archaea` (GTDB) or
  `Eukaryota`/`Viruses` (NCBI), as in sourmash's own combined databases.
- `<prefix>_genome_to_taxid.tsv` gives the NCBI taxid of every genome that has one, the file `yacht convert
  --genome_to_taxid` needs. The GlobDB MAGs have no NCBI taxid.
- `<prefix>_removed_genomes.tsv` lists the genomes `yacht train` left out as equivalent to a genome in the database.
- How to get the sequence of a genome: [retrieving genomes](retrieving_genomes.md).

## Things to know

- **Duplicates are handled by training.** Each source may contain near-identical genomes (e.g. the GenBank and RefSeq
  copies of an assembly in NCBI Virus); `yacht train` keeps one representative per group of genomes within the ANI
  threshold.
- **Viruses at scaled 1000.** In the unions, viruses are sketched at scaled 1000 (downsampled from scaled 50), where a
  30 kb genome has about 30 hashes and many small viruses have none; genomes without any hash are left out. Use the
  scaled-50 viral databases for viruses.
- **Prophages.** `yacht train` treats a genome contained in a larger one at or above the ANI threshold as equivalent
  to it. A phage that is part of a bacterial genome in the database (as a prophage) is therefore represented by that
  genome in the unions.
- **Licenses.** GTDB and GlobDB data are CC BY-SA 4.0, so the databases built from them (including the unions) are
  CC BY-SA 4.0; the ICTV VMR is CC BY 4.0; NCBI places no restrictions on its data. Each database's `README.md` lists
  its sources and what to cite.

## How the databases were built

The scripts are in [`databases/`](../databases) (branch `optimized_databases`):

- `fetch_sources.sh`: downloads the sources (GTDB, GlobDB, and the sketches prepared by the sourmash team for the
  GenBank, ICTV VMR and NCBI viral genomes) and records their md5 sums.
- `sketch_gtdb_r232.sh`, `sketch_globdb_r232.sh`: sketch the GTDB and GlobDB genomes with
  `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000`.
- `assemble_unions.sh`: downsamples the viral sketches to scaled 1000 (`sourmash sig downsample`) and joins the sources
  (`sourmash sig cat`).
- `train_db.sh`: `yacht train --ksize 31 --ani_thresh <ANI>` for each database and threshold.
- `package_db.py`: writes the Zenodo zip (config, manifest, index, lineages, taxids, removed genomes, data card).
