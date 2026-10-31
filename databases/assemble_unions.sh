#!/bin/bash
# Assemble the k31 inputs of the two union databases (scaled 1000):
#   yacht-union-insdc-2026:  GTDB R232 reps + GenBank fungi + protozoa (2025.04) + NCBI RefSeq viruses (2025.01) + ICTV VMR MSL40
#   yacht-union-globdb-2026: GlobDB r232      + GenBank fungi + protozoa (2025.04) + NCBI RefSeq viruses (2025.01) + ICTV VMR MSL40
# Eukaryotes without the assemblies NCBI flags as contaminated (filter_eukaryotes.sh); viruses: the RefSeq genomes
# (extract_ncbi_viruses_refseq.sh) and the ICTV exemplars, downsampled from scaled 50 to 1000 (sketches left empty are
# skipped by yacht train). Also writes the ICTV lineages in the union's rank columns (superkingdom = Viruses).
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/sourmash; W=$YACHT_DB_WORK/unions; E=$YACHT_DB_WORK/eukaryotes; mkdir -p $W
SM="env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash"
# sourmash adds to an existing output zip instead of replacing it: start from scratch
rm -f $W/ncbi-viruses-refseq-2025.01.k31-sc1000.sig.zip $W/ictv-vmr-msl40.k31-sc1000.sig.zip \
  $W/yacht-union-insdc-2026.k31.sig.zip $W/yacht-union-globdb-2026.k31.sig.zip
echo "$(date -Is) downsample the viruses to scaled 1000"
$SM sig downsample --scaled 1000 $YACHT_DB_WORK/ncbi-viruses-refseq/ncbi-viruses-refseq-2025.01.k31-sc50.sig.zip \
  -o $W/ncbi-viruses-refseq-2025.01.k31-sc1000.sig.zip > $W/downsample_refseq.log 2>&1
$SM sig downsample --scaled 1000 $S/2025-ictv-vmr-msl40/vmr_MSL40_v1.dna-k31.sig.zip \
  -o $W/ictv-vmr-msl40.k31-sc1000.sig.zip > $W/downsample_vmr.log 2>&1
$YACHT_ENV/bin/python -I - $S/2025-ictv-vmr-msl40/vmr_MSL40_v1.lineages.csv $W/vmr_MSL40_v1.union-ranks.lineages.csv <<'PY'
import csv
import sys

ranks = ["phylum", "class", "order", "family", "genus", "species"]
with open(sys.argv[1], newline="") as f_in, open(sys.argv[2], "w", newline="") as f_out:
    w = csv.writer(f_out)
    w.writerow(["ident", "superkingdom"] + ranks)
    for row in csv.DictReader(f_in):
        w.writerow([row["ident"], "Viruses"] + [row[r] for r in ranks])
PY
EUK_VIR="$E/genbank-2025.04-fungi-k31.sig.zip $E/genbank-2025.04-protozoa-k31.sig.zip $W/ncbi-viruses-refseq-2025.01.k31-sc1000.sig.zip $W/ictv-vmr-msl40.k31-sc1000.sig.zip"
echo "$(date -Is) union INSDC"
$SM sig cat $YACHT_DB_WORK/gtdb-r232/gtdb-rs232-reps.k31.sig.zip $EUK_VIR -o $W/yacht-union-insdc-2026.k31.sig.zip > $W/cat_insdc.log 2>&1
echo "$(date -Is) union GlobDB"
$SM sig cat $YACHT_DB_WORK/globdb-r232/globdb-r232.k31.sig.zip $EUK_VIR -o $W/yacht-union-globdb-2026.k31.sig.zip > $W/cat_globdb.log 2>&1
for f in $W/cat_insdc.log $W/cat_globdb.log; do tr '\r' '\n' < $f | grep -E "output .* signatures" | tail -1; done
echo "$(date -Is) ASSEMBLE_DONE"
