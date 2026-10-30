# Building the YACHT 2.0 pretrained databases

These scripts build the databases described in [`docs/databases.md`](../docs/databases.md) and stage them for the
YACHT Zenodo community. They use only existing tools and identifiers: the source genomes or the sketches the sourmash
team prepared, `sourmash scripts manysketch` (sourmash_plugin_branchwater) to sketch FASTA files, `sourmash sig
downsample/cat/extract`, and `yacht train`.

## Setup

```bash
conda activate yacht                  # an environment with YACHT >= 2.0 (env/yacht_env.yml + pip install .)
export YACHT_DB_WORK=/big/disk/dbbuild  # several TB of free space (see below)
```

`config.sh` reads `YACHT_DB_WORK` (default `./dbbuild`) and `YACHT_ENV` (default: the active conda environment).

## Steps

```bash
bash databases/fetch_sources.sh sourmash   # sourmash-prepared sketches: GTDB rs226, GenBank 2025.04 fungi/protozoa, ICTV VMR MSL40, NCBI viruses 2025.01
bash databases/fetch_sources.sh gtdb       # GTDB R232: species representative genomes (192 GB) + metadata
bash databases/fetch_sources.sh globdb     # GlobDB r232: genomes (292 GB), taxonomy, ID dictionaries
bash databases/fetch_sources.sh ncbi       # NCBI GenBank assembly summaries of fungi and protozoa (contamination flags)
bash databases/filter_eukaryotes.sh        # fungi/protozoa sketches without the assemblies NCBI flags as contaminated
bash databases/extract_ncbi_viruses_refseq.sh   # the RefSeq (GCF_) genomes of the NCBI viral sketches
bash databases/sketch_gtdb_r232.sh         # md5 check, extract, manysketch k=21,31,51 (names: "<accession> <NCBI organism name>")
bash databases/sketch_globdb_r232.sh       # md5 check, extract the non-GTDB genomes, manysketch, + GTDB R232 sketches
                                           # (both call rezip_manysketch.sh: sourmash 4.9 cannot open manysketch zips > 4 GiB)
bash databases/assemble_unions.sh          # RefSeq + ICTV viral sketches downsampled to scaled 1000; sourmash sig cat of the parts
# training: one line "REF_ZIP PREFIX ANI" per database and ANI threshold, e.g.
#   $YACHT_DB_WORK/gtdb-r232/gtdb-rs232-reps.k31.sig.zip gtdb-rs232-reps.k31_0.95 0.95
bash databases/run_jobs.sh jobs.txt 3 48   # at most 3 trainings at a time, 48 threads each
bash databases/package_all.sh              # zips with config, manifest, index, lineages, taxids, removed genomes, data card
bash databases/verify_db.sh $YACHT_DB_WORK/package/<prefix>_pretrained.zip sample1.sig.zip ...   # unzip as a user, yacht run
python databases/stage_zenodo.py $YACHT_DB_WORK/package $YACHT_DB_WORK/zenodo   # metadata + upload/publish scripts
python databases/collect_cards.py $YACHT_DB_WORK/package docs/database_cards     # the data cards, for the repository
```

### Uploading to Zenodo

`stage_zenodo.py` writes one folder per database (zip, data card `README.md`, `metadata.json`, `md5.txt`), `INDEX.tsv`
and two scripts, to be run by an owner of the [YACHT community](https://zenodo.org/communities/yacht):

1. Create a personal access token on Zenodo (account menu > Applications > Personal access tokens > New token,
   https://zenodo.org/account/settings/applications/tokens/new/) with the scope `deposit:write`, and `deposit:actions`
   to publish with the script.
2. `ZENODO_TOKEN=<token> bash upload_to_zenodo.sh` creates one **draft** per database, with the metadata (title,
   license, related DOIs, keywords) and the request to include it in the `yacht` community, uploads the zip and the
   data card, and checks the zip's md5 on Zenodo. Re-running it skips the databases that already have a draft.
   With `ZENODO_API=https://sandbox.zenodo.org/api` (and a sandbox token and `ZENODO_COMMUNITY=<a sandbox
   community>`) everything can be tried on the Zenodo sandbox first.
3. Review the drafts on https://zenodo.org/me/uploads, then publish them there or with
   `ZENODO_TOKEN=<token> bash publish_zenodo_drafts.sh` (publishing mints the DOIs and cannot be undone).
4. A record that asks to join a community is included once a community curator accepts the request (the community's
   "Requests" tab). Then `yacht download pretrained_ref_db --list` shows the database.

`databases/db_sources.py` holds the data-card texts (sources, licenses, citations, identifier schemes);
`databases/remote_manifest.py` summarizes the manifest of a remote sourmash zip without downloading it.

## Cost (measured on a 768-core / 3 TB server, October 2026)

Steps before training (wall time):

| step | time | notes |
|---|---:|---|
| download GTDB R232 genomes (192 GB) / GlobDB r232 genomes (292 GB) | ~3.5 h each, in parallel | 17-23 MB/s |
| sketch GTDB R232 (199,923 genomes; manysketch, 64 threads) | 28 min | + 33 min `rezip_manysketch.sh` |
| sketch the GlobDB MAGs (146,310 genomes; 64 threads) | 14 min | + 18 min `rezip_manysketch.sh` |
| `filter_eukaryotes.sh` | 14 min | |
| assemble the two unions (`sourmash sig cat`) | 25 + 30 min | |

Training (`yacht train --num_threads 48`, up to 4 at a time) and the resulting databases (104 GB of zips in total):

| database | ANI | sketches trained on | genomes in the database | zip (GB) | `yacht train` wall time | peak RAM (GiB) |
|---|---|---:|---:|---:|---:|---:|
| `genbank-2025.04-fungi-k31` | 0.80 | 21,070 | 967 | 0.35 | 6:13.65 | 26 |
| `genbank-2025.04-fungi-k31` | 0.95 | 21,070 | 4,990 | 1.35 | 6:29.53 | 26 |
| `genbank-2025.04-fungi-k31` | 0.995 | 21,070 | 9,081 | 1.90 | 6:45.20 | 26 |
| `genbank-2025.04-fungi-k31` | 0.9995 | 21,070 | 15,564 | 2.30 | 6:52.99 | 26 |
| `genbank-2025.04-protozoa-k31` | 0.80 | 2,346 | 237 | 0.26 | 1:40.87 | 8 |
| `genbank-2025.04-protozoa-k31` | 0.95 | 2,346 | 891 | 0.54 | 1:43.06 | 8 |
| `genbank-2025.04-protozoa-k31` | 0.995 | 2,346 | 1,512 | 0.61 | 1:45.66 | 8 |
| `genbank-2025.04-protozoa-k31` | 0.9995 | 2,346 | 2,076 | 0.63 | 1:49.34 | 8 |
| `globdb-r232.k31` | 0.80 | 346,233 | 17,558 | 0.59 | 35:40.52 | 117 |
| `globdb-r232.k31` | 0.95 | 346,233 | 298,599 | 7.11 | 34:58.13 | 117 |
| `globdb-r232.k31` | 0.995 | 346,233 | 345,661 | 7.79 | 39:07.26 | 117 |
| `globdb-r232.k31` | 0.9995 | 346,233 | 346,188 | 7.79 | 39:55.78 | 117 |
| `gtdb-rs232-reps.k31` | 0.80 | 199,923 | 14,568 | 0.49 | 17:30.18 | 84 |
| `gtdb-rs232-reps.k31` | 0.95 | 199,923 | 184,616 | 4.83 | 18:06.56 | 84 |
| `gtdb-rs232-reps.k31` | 0.995 | 199,923 | 199,923 | 5.14 | 22:26.36 | 84 |
| `gtdb-rs232-reps.k31` | 0.9995 | 199,923 | 199,923 | 5.14 | 22:06.69 | 84 |
| `ictv-vmr-msl40.k31-sc50` | 0.80 | 17,764 | 8,987 | 0.06 | 0:12.91 | 2 |
| `ictv-vmr-msl40.k31-sc50` | 0.95 | 17,764 | 13,530 | 0.11 | 0:14.86 | 2 |
| `ictv-vmr-msl40.k31-sc50` | 0.995 | 17,764 | 17,530 | 0.14 | 0:13.21 | 2 |
| `ictv-vmr-msl40.k31-sc50` | 0.9995 | 17,764 | 17,694 | 0.15 | 0:13.06 | 2 |
| `ncbi-viruses-2025.01.k31-sc50` | 0.80 | 230,967 | 19,218 | 0.29 | 5:41.74 | 72 |
| `ncbi-viruses-2025.01.k31-sc50` | 0.95 | 230,967 | 35,447 | 0.48 | 2:49.02 | 42 |
| `ncbi-viruses-2025.01.k31-sc50` | 0.995 | 230,967 | 73,235 | 0.66 | 2:04.43 | 34 |
| `ncbi-viruses-2025.01.k31-sc50` | 0.9995 | 230,967 | 127,987 | 0.74 | 2:01.33 | 33 |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | 0.80 | 16,279 | 8,255 | 0.06 | 0:12.44 | 2 |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | 0.95 | 16,279 | 12,227 | 0.11 | 0:12.71 | 2 |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | 0.995 | 16,279 | 15,797 | 0.13 | 0:13.07 | 2 |
| `ncbi-viruses-refseq-2025.01.k31-sc50` | 0.9995 | 16,279 | 16,081 | 0.14 | 0:13.44 | 2 |
| `yacht-union-globdb-2026.k31` | 0.80 | 402,672 | 31,092 | 1.01 | 41:34.17 | 145 |
| `yacht-union-globdb-2026.k31` | 0.95 | 402,672 | 319,149 | 9.06 | 42:11.31 | 145 |
| `yacht-union-globdb-2026.k31` | 0.995 | 402,672 | 376,021 | 10.35 | 42:45.08 | 145 |
| `yacht-union-globdb-2026.k31` | 0.9995 | 402,672 | 383,918 | 10.77 | 43:30.05 | 145 |
| `yacht-union-insdc-2026.k31` | 0.80 | 256,362 | 28,163 | 0.91 | 31:11.52 | 113 |
| `yacht-union-insdc-2026.k31` | 0.95 | 256,362 | 205,219 | 6.65 | 30:55.26 | 113 |
| `yacht-union-insdc-2026.k31` | 0.995 | 256,362 | 230,291 | 7.56 | 33:46.23 | 113 |
| `yacht-union-insdc-2026.k31` | 0.9995 | 256,362 | 237,662 | 8.12 | 31:13.23 | 113 |

`yacht train`'s memory is dominated by the C++ hash index of the comparison core (about 140 bytes per reference
hash), since its comparison matrix is limited to 32 GiB. The intermediate files of a training (up to about 70 GB for
the unions) can be deleted after packaging. Peak disk use of the whole build was about 1.3 TB (downloads 463 GB,
extracted genomes about 280 GB, training intermediates about 500 GB).
