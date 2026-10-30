#!/bin/bash
# Sketch the GlobDB r232 genomes that are not GTDB R232 species representatives (the MAGs of the other catalogues;
# k=21,31,51; scaled=1000), and assemble the GlobDB r232 k31 database input: the GTDB R232 sketches (same genomes,
# versioned accessions) + the MAG sketches. Can run while sketch_gtdb_r232.sh runs; its last step waits for it.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/globdb-r232; W=$YACHT_DB_WORK/globdb-r232; G=$YACHT_DB_WORK/gtdb-r232; mkdir -p $W
THREADS=${THREADS:-64}
YACHT_CHECK_PY="$YACHT_ENV/bin/python -I"
echo "$(date -Is) md5 check"
want=$(awk '$2 ~ /globdb_r232_genome_fasta.tar.gz$/ {print $1}' $S/globdb_r232_md5sum.txt)
got=$(md5sum $S/globdb_r232_genome_fasta.tar.gz | cut -d' ' -f1)
echo "md5 want=$want got=$got"; [ "$want" = "$got" ]
echo "$(date -Is) extract (genomes without a GCA/GCF accession; in parallel, three GTDB genomes to check that GlobDB's copies are the same)"
CHECK="GCF_000020485 GCF_000020525 GCF_000020985"
mkdir -p $W/check
( tar xzf $S/globdb_r232_genome_fasta.tar.gz -C $W/check --wildcards $(for c in $CHECK; do echo "*/$c.fa.gz"; done) ) &
check_pid=$!
[ -f $W/.extracted ] || { tar xzf $S/globdb_r232_genome_fasta.tar.gz -C $W --exclude '*/GCA_*' --exclude '*/GCF_*' && touch $W/.extracted; }
wait $check_pid
for c in $CHECK; do
  f=$(find $W/check -name "$c.fa.gz"); g=$(find $G/gtdb_genomes_reps_r232 -name "${c}.*_genomic.fna.gz")
  $YACHT_CHECK_PY - "$f" "$g" <<'PY' || echo "check of $c failed"
import sys
import screed
from sourmash import MinHash
def sketch(path):
    mh = MinHash(0, 31, scaled=1000)
    for record in screed.open(path):
        mh.add_sequence(record.sequence, force=True)
    return mh
a, b = sketch(sys.argv[1]), sketch(sys.argv[2])
same = "IDENTICAL sketch" if a.hashes.keys() == b.hashes.keys() else "DIFFERENT"
print(f"GlobDB {sys.argv[1].split('/')[-1]} vs GTDB {sys.argv[2].split('/')[-1]}: {same} ({len(a)} hashes)")
PY
done
echo "$(date -Is) manysketch CSV, lineages; GTDB part check"
$YACHT_ENV/bin/python -I $SCRIPTS/globdb_sketch_csv.py $S $W/globdb_r232_genome_fasta $G/gtdb-rs232-reps.lineages.csv $W/globdb-r232-mags.manysketch.csv
echo "$(date -Is) manysketch"
[ -f $W/.manysketch_done ] || { env PATH=$YACHT_ENV/bin:$PATH RAYON_NUM_THREADS=$THREADS $YACHT_ENV/bin/sourmash scripts manysketch \
  $W/globdb-r232-mags.manysketch.csv -p dna,k=21,k=31,k=51,scaled=1000 -c $THREADS -o $W/globdb-r232-mags.manysketch.zip \
  > $W/manysketch.log 2>&1 && touch $W/.manysketch_done; }
echo "$(date -Is) rewrite with sourmash sig cat (all k, k31)"
bash $SCRIPTS/rezip_manysketch.sh $W/globdb-r232-mags.manysketch.zip $W/globdb-r232-mags.k21-k31-k51.sig.zip \
  $W/globdb-r232-mags.k31.sig.zip
rm $W/globdb-r232-mags.manysketch.zip
echo "$(date -Is) k31: GTDB R232 + MAGs (waiting for sketch_gtdb_r232.sh)"
until [ -f $G/.sketch_done ]; do sleep 60; done
rm -f $W/globdb-r232.k31.sig.zip  # sourmash adds to an existing output zip instead of replacing it
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig cat $G/gtdb-rs232-reps.k31.sig.zip $W/globdb-r232-mags.k31.sig.zip -o $W/globdb-r232.k31.sig.zip > $W/cat_k31.log 2>&1
echo "$(date -Is) SKETCH_DONE"
