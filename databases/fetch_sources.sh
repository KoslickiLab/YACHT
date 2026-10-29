#!/bin/bash
# Fetch database sources (resumable). Usage: fetch_sources.sh {gtdb|globdb|sourmash}
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
S=$YACHT_DB_WORK/src
get() { echo "$(date -Is) START $1"; curl -sS -L -C - --retry 20 --retry-delay 30 --retry-all-errors -o "$2" "$1"; echo "$(date -Is) END rc=$? $(stat -c %s "$2" 2>/dev/null) $2"; }
case $1 in
  gtdb)
    G=https://data.gtdb.ecogenomic.org/releases/release232/232.0
    for f in VERSION.txt RELEASE_NOTES.txt FILE_DESCRIPTIONS.txt MD5SUM.txt METHODS.txt bac120_metadata_r232.tsv.gz ar53_metadata_r232.tsv.gz bac120_taxonomy_r232.tsv.gz ar53_taxonomy_r232.tsv.gz; do get $G/$f $S/gtdb-r232/$f; done
    get $G/genomic_files_reps/gtdb_genomes_reps_r232.tar.gz $S/gtdb-r232/gtdb_genomes_reps_r232.tar.gz;;
  globdb)
    G=https://fileshare.lisc.univie.ac.at/globdb/globdb_r232
    for f in globdb_r232_md5sum.txt globdb_r232_taxonomy.tsv.gz globdb_r232_dictionaries.tar.gz globdb_r232_dataset_list.tsv; do get $G/$f $S/globdb-r232/$f; done
    get $G/globdb_r232_genome_fasta.tar.gz $S/globdb-r232/globdb_r232_genome_fasta.tar.gz;;
  sourmash)
    B=https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db
    for f in gtdb-rs226/gtdb-rs226-reps.k31.sig.zip gtdb-rs226/gtdb-rs226-reps.lineages.csv \
             genbank-2025.04/genbank-20250408-fungi-k31.zip genbank-2025.04/genbank-20250408-protozoa-k31.zip \
             genbank-2025.04/lineages.fungi.csv genbank-2025.04/lineages.protozoa.csv \
             2025-ictv-vmr-msl40/vmr_MSL40_v1.dna-k31.sig.zip 2025-ictv-vmr-msl40/vmr_MSL40_v1.lineages.csv \
             ncbi-viruses-2025.01/ncbi-viruses.dna.k=31.scaled=50.sig.zip ncbi-viruses-2025.01/ncbi-viruses.lineages.csv; do
      mkdir -p $S/sourmash/$(dirname $f); get "$B/$f" "$S/sourmash/$f"; done;;
esac
echo "$(date -Is) FETCH_DONE $1"
