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

**These databases need YACHT 2.0.0 or later.** They hold the YACHT 2.0 index but not the genome sketches that YACHT
1.x reads. This is recorded in each database's `<prefix>_config.json` (`"min_yacht_version": "2.0.0"`, with the YACHT
version that made it as `yacht_version`), so YACHT refuses a database that needs a newer version than the one
installed, and it is stated in each data card and Zenodo record (title: `... (for YACHT >= 2.0.0)`). The data cards of
all databases are in [database_cards](database_cards/README.md).

| name | content | genomes before training | identifiers | taxonomy | scaled |
|---|---|---:|---|---|---:|
| `gtdb-rs232-reps.k31` | GTDB R232 species representatives (bacteria, archaea) | 199,923 | versioned GCA/GCF accessions | GTDB R232 | 1000 |
| `globdb-r232.k31` | GlobDB r232: GTDB R232 species representatives + species-level MAGs of 25 catalogues | 346,233 | GCA/GCF accessions + GlobDB IDs | GlobDB r232 (GTDB R232) | 1000 |
| `genbank-2025.04-fungi-k31` | GenBank fungal genomes, April 2025, without the assemblies NCBI flags as contaminated | 21,070 | GCA accessions | NCBI | 1000 |
| `genbank-2025.04-protozoa-k31` | GenBank protozoan genomes, April 2025, without the assemblies NCBI flags as contaminated | 2,346 | GCA accessions | NCBI | 1000 |
| `ictv-vmr-msl40.k31-sc50` | ICTV exemplar viruses (VMR MSL40) | 17,765 | GCA accessions | ICTV MSL40 | 50 |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | NCBI RefSeq viral genomes, January 2025 (curated) | 16,280 | GCF accessions | NCBI | 50 |
| `ncbi-viruses-2025.01.k31-sc50` | all NCBI viral genome assemblies (GenBank + RefSeq), January 2025; see the note below | 230,973 | GCA/GCF accessions | NCBI | 50 |
| `yacht-union-insdc-2026.k31` | GTDB R232 reps + GenBank fungi + protozoa (both filtered as above) + NCBI RefSeq viruses + ICTV VMR exemplars | 257,384 | GCA/GCF accessions only | GTDB (bacteria, archaea), NCBI (eukaryotes, RefSeq viruses), ICTV (VMR exemplars) | 1000 |
| `yacht-union-globdb-2026.k31` | GlobDB r232 + GenBank fungi + protozoa (both filtered as above) + NCBI RefSeq viruses + ICTV VMR exemplars | 403,694 | GCA/GCF accessions + GlobDB IDs | GlobDB (bacteria, archaea), NCBI (eukaryotes, RefSeq viruses), ICTV (VMR exemplars) | 1000 |

The numbers of genomes in each database after training are in the database's `README.md`; [`databases/README.md`](../databases/README.md) lists them for all databases and thresholds, with the zip sizes and build costs.

## Which database?

- **Bacteria and archaea, easiest to follow up:** `gtdb-rs232-reps.k31`. Every genome has an NCBI accession.
- **Bacteria and archaea, most comprehensive:** `globdb-r232.k31`. It adds about 146,000 species that are only known
  from metagenome-assembled genomes. Their genomes are retrieved from the GlobDB release instead of NCBI (see
  [retrieving genomes](retrieving_genomes.md)).
- **Everything (bacteria, archaea, fungi, protists, viruses):** `yacht-union-insdc-2026.k31` (every genome retrievable
  from NCBI by its accession; the easiest to follow up) or `yacht-union-globdb-2026.k31` (most comprehensive; GlobDB
  genomes are retrieved from the GlobDB release).
- **Viruses:** the scaled-50 databases `ncbi-viruses-refseq-2025.01.k31-sc50` (curated NCBI RefSeq genomes) and
  `ictv-vmr-msl40.k31-sc50` (curated exemplars with ICTV taxonomy). Sketch the sample with `--scaled 50`; the same
  sketch also works with the scaled-1000 databases, since YACHT downsamples it. `ncbi-viruses-2025.01.k31-sc50` holds
  all NCBI viral assemblies, including GenBank assemblies that are not single viral genomes (see below).

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

- **ANI 0.80 is very coarse, much coarser than genus level.** At k = 31 it corresponds to a k-mer containment of
  0.80^31 = 0.001, about 4 shared hashes for a 4 Mbp genome at scaled 1000. Overlaps this small also occur between
  distant genomes: in `gtdb-rs232-reps.k31_0.80`, 70% of the GTDB genera, 41% of the families and 15 phyla keep no
  genome of their own, and 13% of the left-out genomes that have a genome of the database at ANI >= 0.80 are
  represented by a genome of another phylum (at ANI 0.95: every family keeps a genome, and 99.9% of left-out genomes
  are represented by a genome of the same genus). Each database's `README.md` has this table ("How coarse is this
  database?"). Use 0.95 or higher for genus- or species-level results.
- **Equivalence chains.** `yacht train` goes from the smallest to the largest genome and leaves a genome out when it is
  equivalent to a larger genome that has not been left out; that genome can be left out later in turn. A genome left
  out through such a chain can have no genome of the database at the ANI threshold (an empty `closest_genome_in_database`
  in `<prefix>_removed_genomes.tsv`).
- **Assemblies that are not what their label says.** `yacht train` keeps the larger of two equivalent genomes, so an
  assembly that contains other genomes can represent them. In GenBank, some "viral" assemblies are collections of many
  viral genomes (`GCA_024181985.1` "uncultured virus", 334 Mbp) or contain their host's genome (`GCA_038430685.1`
  "Enterococcus phage Efs 1-1-1", 3 Mbp), and some eukaryotic assemblies contain bacterial sequence. The unions
  therefore use only curated viral genomes (RefSeq, ICTV exemplars) and leave out the fungal and protozoan assemblies
  that NCBI flags as contaminated; `ncbi-viruses-2025.01.k31-sc50` keeps all viral assemblies and says so in its
  `README.md`. Other contaminated eukaryotic assemblies may remain, and a few ICTV exemplars are sketches of the host
  assembly that carries the virus (prophages, gene transfer agents, e.g. `GCA_000091085.2` "Chlamydia phage phiCPAR39",
  the 1.25 Mbp *Chlamydia pneumoniae* AR39 genome); `<prefix>_removed_genomes.tsv` shows which genome represents each
  left-out genome.
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
  `sourmash scripts manysketch -p dna,k=21,k=31,k=51,scaled=1000`. sourmash 4.9 cannot open the zips over 4 GiB that
  manysketch writes, so `rezip_manysketch.sh` unpacks them and writes them again with `sourmash sig cat` (all
  k-sizes, and k31 only for training), checking that the sketches are the same.
- `filter_eukaryotes.sh`: leaves out the GenBank fungal and protozoan assemblies that NCBI flags as contaminated
  (`excluded_from_refseq` in the NCBI assembly summaries; `sourmash sig extract --picklist ...:exclude`).
- `extract_ncbi_viruses_refseq.sh`: selects the RefSeq (`GCF_`) genomes of the NCBI viral sketches
  (`sourmash sig extract --picklist`).
- `assemble_unions.sh`: downsamples the RefSeq and ICTV viral sketches to scaled 1000 (`sourmash sig downsample`) and
  joins the sources (`sourmash sig cat`).
- `train_db.sh`: `yacht train --ksize 31 --ani_thresh <ANI>` for each database and threshold.
- `package_db.py`: writes the Zenodo zip (config, manifest, index, lineages, taxids, removed genomes, data card).
