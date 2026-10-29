#!/bin/bash
# Sketch the GTDB R232 species representative genomes (k=21,31,51; scaled=1000; no abundances, like sourmash's
# own GTDB databases), named "<versioned accession> <NCBI organism name>" as in sourmash's GTDB databases.
# Input: $YACHT_DB_WORK/src/gtdb-r232/gtdb_genomes_reps_r232.tar.gz (+ bac120/ar53 metadata), checked against MD5SUM.txt.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/gtdb-r232; W=$YACHT_DB_WORK/gtdb-r232; mkdir -p $W
THREADS=${THREADS:-64}
echo "$(date -Is) md5 check"
want=$(awk '$2 ~ /gtdb_genomes_reps_r232.tar.gz$/ {print $1}' $S/MD5SUM.txt)
got=$(md5sum $S/gtdb_genomes_reps_r232.tar.gz | cut -d' ' -f1)
echo "md5 want=$want got=$got"; [ "$want" = "$got" ]
echo "$(date -Is) extract"
[ -f $W/.extracted ] || { tar xzf $S/gtdb_genomes_reps_r232.tar.gz -C $W && touch $W/.extracted; }
echo "$(date -Is) manysketch CSV"
$YACHT_ENV/bin/python -I $SCRIPTS/gtdb_sketch_csv.py $S $W/gtdb_genomes_reps_r232 $W/gtdb-rs232-reps.manysketch.csv
echo "$(date -Is) manysketch"
env PATH=$YACHT_ENV/bin:$PATH RAYON_NUM_THREADS=$THREADS $YACHT_ENV/bin/sourmash scripts manysketch $W/gtdb-rs232-reps.manysketch.csv \
  -p dna,k=21,k=31,k=51,scaled=1000 -c $THREADS -o $W/gtdb-rs232-reps.k21-k31-k51.sig.zip > $W/manysketch.log 2>&1
echo "$(date -Is) extract k31"
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig extract --ksize 31 $W/gtdb-rs232-reps.k21-k31-k51.sig.zip -o $W/gtdb-rs232-reps.k31.sig.zip > $W/extract_k31.log 2>&1
echo "$(date -Is) SKETCH_DONE"
