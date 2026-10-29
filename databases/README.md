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
bash databases/sketch_gtdb_r232.sh         # md5 check, extract, manysketch k=21,31,51 (names: "<accession> <NCBI organism name>")
bash databases/sketch_globdb_r232.sh       # md5 check, extract the non-GTDB genomes, manysketch, + GTDB R232 sketches
bash databases/assemble_unions.sh          # viral sketches downsampled to scaled 1000; sourmash sig cat of the parts
# training: one line "REF_ZIP PREFIX ANI" per database and ANI threshold, e.g.
#   $YACHT_DB_WORK/gtdb-r232/gtdb-rs232-reps.k31.sig.zip gtdb-rs232-reps.k31_0.95 0.95
bash databases/run_jobs.sh jobs.txt 3 48   # at most 3 trainings at a time, 48 threads each
bash databases/package_all.sh              # zips with config, manifest, index, lineages, taxids, removed genomes, data card
python databases/stage_zenodo.py $YACHT_DB_WORK/package $YACHT_DB_WORK/zenodo   # metadata + upload script (drafts only)
```

`databases/db_sources.py` holds the data-card texts (sources, licenses, citations, identifier schemes);
`databases/remote_manifest.py` summarizes the manifest of a remote sourmash zip without downloading it.

## Cost (measured on a 768-core / 3 TB server, 48-64 threads per job)

See the table in the pull request / build notes; `yacht train`'s memory is dominated by the C++ hash index of the
comparison core (about 140 bytes per reference hash), since its comparison matrix is limited to 32 GiB.
