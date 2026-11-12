#!/bin/bash
# The NCBI RefSeq viral genomes (GCF_ accessions) of the sourmash team's NCBI viruses 2025.01 sketches (k31, scaled 50):
# the curated subset, without the GenBank-only multi-genome collections and host-containing assemblies.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src/sourmash/ncbi-viruses-2025.01; W=$YACHT_DB_WORK/ncbi-viruses-refseq; mkdir -p $W
$YACHT_ENV/bin/python -I - $S/ncbi-viruses.dna.k=31.scaled=50.sig.zip $W/refseq.picklist.csv <<'PY'
import csv
import sys
import zipfile

with zipfile.ZipFile(sys.argv[1]) as z:
    lines = z.read("SOURMASH-MANIFEST.csv").decode().splitlines()
idents = sorted({r["name"].split()[0] for r in csv.DictReader(line for line in lines if not line.startswith("#"))
                 if r["name"].startswith("GCF_")})
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ident"])
    w.writerows([i] for i in idents)
print(f"RefSeq (GCF_) viral assemblies: {len(idents)}")
PY
rm -f $W/ncbi-viruses-refseq-2025.01.k31-sc50.sig.zip  # sourmash adds to an existing output zip instead of replacing it
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig extract $S/ncbi-viruses.dna.k=31.scaled=50.sig.zip \
  --picklist $W/refseq.picklist.csv:ident:ident -o $W/ncbi-viruses-refseq-2025.01.k31-sc50.sig.zip > $W/extract.log 2>&1
tr '\r' '\n' < $W/extract.log | grep -E "extracted|picklist" | tail -3
echo "$(date -Is) EXTRACT_DONE"
