#!/bin/bash
# Leave out the GenBank fungal and protozoan assemblies that NCBI flags as contaminated (`excluded_from_refseq`
# contains "contaminated" in the current or historical NCBI assembly summary; fetch_sources.sh ncbi) from the
# sourmash team's GenBank 2025.04 k31 sketches, with `sourmash sig extract --picklist ...:exclude`.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
S=$YACHT_DB_WORK/src; W=$YACHT_DB_WORK/eukaryotes; mkdir -p $W
for org in fungi protozoa; do
  $YACHT_ENV/bin/python -I - $S/ncbi-assembly-summaries $org $S/sourmash/genbank-2025.04/genbank-20250408-$org-k31.zip \
    $W/$org.ncbi-contaminated.csv <<'PY'
import csv
import sys
import zipfile

summaries, org, sketches, out = sys.argv[1:5]
flags = {}
for kind in ("assembly_summary_historical", "assembly_summary"):  # current flags win
    with open(f"{summaries}/genbank_{org}_{kind}.txt") as f:
        next(f)
        header = next(f).lstrip("#").rstrip("\n").split("\t")
        for line in f:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            flags[row["assembly_accession"]] = row["excluded_from_refseq"]
with zipfile.ZipFile(sketches) as z:
    lines = z.read("SOURMASH-MANIFEST.csv").decode().splitlines()
names = {r["name"].split()[0]: r["name"] for r in csv.DictReader(line for line in lines if not line.startswith("#"))}
flagged = sorted(i for i in names if "contaminated" in flags.get(i, "").lower())
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ident", "name", "excluded_from_refseq"])
    w.writerows([i, names[i], flags[i]] for i in flagged)
print(f"{org}: {len(names)} sketches, {sum(1 for i in names if i not in flags)} not in the NCBI summaries, "
      f"{len(flagged)} flagged as contaminated")
PY
  rm -f $W/genbank-2025.04-$org-k31.sig.zip  # sourmash adds to an existing output zip instead of replacing it
  env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig extract $S/sourmash/genbank-2025.04/genbank-20250408-$org-k31.zip \
    --picklist $W/$org.ncbi-contaminated.csv:ident:ident:exclude -o $W/genbank-2025.04-$org-k31.sig.zip > $W/extract_$org.log 2>&1
  tr '\r' '\n' < $W/extract_$org.log | grep -E "extracted" | tail -1
done
echo "$(date -Is) FILTER_DONE"
