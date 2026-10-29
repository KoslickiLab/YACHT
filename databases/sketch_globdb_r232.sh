#!/bin/bash
# Sketch the GlobDB r232 genomes that are not GTDB R232 species representatives (the MAGs of the other catalogues;
# k=21,31,51; scaled=1000), and assemble the GlobDB r232 k31 database input: the GTDB R232 sketches (same genomes,
# versioned accessions) + the MAG sketches. Needs sketch_gtdb_r232.sh to have run.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/globdb-r232; W=$YACHT_DB_WORK/globdb-r232; G=$YACHT_DB_WORK/gtdb-r232; mkdir -p $W
THREADS=${THREADS:-64}
echo "$(date -Is) md5 check"
want=$(awk '$2 ~ /globdb_r232_genome_fasta.tar.gz$/ {print $1}' $S/globdb_r232_md5sum.txt)
got=$(md5sum $S/globdb_r232_genome_fasta.tar.gz | cut -d' ' -f1)
echo "md5 want=$want got=$got"; [ "$want" = "$got" ]
echo "$(date -Is) extract (genomes without a GCA/GCF accession)"
[ -f $W/.extracted ] || { tar xzf $S/globdb_r232_genome_fasta.tar.gz -C $W --exclude '*/GCA_*' --exclude '*/GCF_*' && touch $W/.extracted; }
echo "$(date -Is) manysketch CSV, lineages; GTDB part check"
$YACHT_ENV/bin/python -I $SCRIPTS/globdb_sketch_csv.py $S $W/globdb_r232_genome_fasta $G/gtdb-rs232-reps.lineages.csv $W/globdb-r232-mags.manysketch.csv
echo "$(date -Is) manysketch"
env PATH=$YACHT_ENV/bin:$PATH RAYON_NUM_THREADS=$THREADS $YACHT_ENV/bin/sourmash scripts manysketch $W/globdb-r232-mags.manysketch.csv \
  -p dna,k=21,k=31,k=51,scaled=1000 -c $THREADS -o $W/globdb-r232-mags.k21-k31-k51.sig.zip > $W/manysketch.log 2>&1
echo "$(date -Is) k31: GTDB R232 + MAGs"
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig extract --ksize 31 $W/globdb-r232-mags.k21-k31-k51.sig.zip -o $W/globdb-r232-mags.k31.sig.zip > $W/extract_k31.log 2>&1
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig cat $G/gtdb-rs232-reps.k31.sig.zip $W/globdb-r232-mags.k31.sig.zip -o $W/globdb-r232.k31.sig.zip > $W/cat_k31.log 2>&1
echo "$(date -Is) SKETCH_DONE"
