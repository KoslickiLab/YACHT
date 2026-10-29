#!/bin/bash
# Assemble the k31 inputs of the two union databases (scaled 1000):
#   yacht-union-insdc-2026:  GTDB R232 reps + GenBank fungi + protozoa (2025.04) + NCBI viruses (2025.01)
#   yacht-union-globdb-2026: GlobDB r232      + GenBank fungi + protozoa (2025.04) + NCBI viruses (2025.01)
# The viral sketches (scaled 50) are downsampled to scaled 1000 first; sketches left empty are skipped by yacht train.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/sourmash; W=$YACHT_DB_WORK/unions; mkdir -p $W
SM="env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash"
echo "$(date -Is) downsample NCBI viruses to scaled 1000"
[ -f $W/ncbi-viruses-2025.01.k31-sc1000.sig.zip ] || $SM sig downsample --scaled 1000 \
  $S/ncbi-viruses-2025.01/ncbi-viruses.dna.k=31.scaled=50.sig.zip -o $W/ncbi-viruses-2025.01.k31-sc1000.sig.zip > $W/downsample.log 2>&1
EUK_VIR="$S/genbank-2025.04/genbank-20250408-fungi-k31.zip $S/genbank-2025.04/genbank-20250408-protozoa-k31.zip $W/ncbi-viruses-2025.01.k31-sc1000.sig.zip"
echo "$(date -Is) union INSDC"
$SM sig cat $YACHT_DB_WORK/gtdb-r232/gtdb-rs232-reps.k31.sig.zip $EUK_VIR -o $W/yacht-union-insdc-2026.k31.sig.zip > $W/cat_insdc.log 2>&1
echo "$(date -Is) union GlobDB"
$SM sig cat $YACHT_DB_WORK/globdb-r232/globdb-r232.k31.sig.zip $EUK_VIR -o $W/yacht-union-globdb-2026.k31.sig.zip > $W/cat_globdb.log 2>&1
echo "$(date -Is) ASSEMBLE_DONE"
