#!/bin/bash
# sourmash 4.9 cannot open the zip files over 4 GiB that `sourmash scripts manysketch` writes ("Invalid Zip archive:
# Extra data field contains disk number"; Python's zipfile and unzip read them), while it reads the large zips it
# writes itself. Unpack the manysketch zip and let `sourmash sig cat` write the collection again, once with all
# k-sizes and once with k31 only (the training input), then check the counts against the manysketch manifest.
# Usage: rezip_manysketch.sh MANYSKETCH_ZIP OUT_ALL_K_ZIP OUT_K31_ZIP
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
in=$1; all=$2; k31=$3
d=$(dirname $all)/$(basename $in .zip).unpacked
rm -rf "${d:?}"; mkdir -p $d
rm -f $all $k31  # sourmash adds to an existing output zip instead of replacing it
echo "$(date -Is) unzip $in"
unzip -q $in 'signatures/*' -d $d
echo "$(date -Is) sourmash sig cat: $(ls $d/signatures | wc -l) sketch files"
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig cat -k 31 $d -o $k31 > $k31.log 2>&1 &
k31_pid=$!
env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/sourmash sig cat $d -o $all > $all.log 2>&1
wait $k31_pid
$YACHT_ENV/bin/python -I - "$in" "$all" "$k31" <<'PY'
import collections
import csv
import sys
import zipfile


def manifest(path):
    with zipfile.ZipFile(path) as z:
        lines = z.read("SOURMASH-MANIFEST.csv").decode().splitlines()
    rows = list(csv.DictReader(line for line in lines if not line.startswith("#")))
    return collections.Counter((r["ksize"], r["name"], r["md5"]) for r in rows)


source, every_k, k31 = (manifest(p) for p in sys.argv[1:4])
want_k31 = collections.Counter({key: n for key, n in source.items() if key[0] == "31"})
print(f"sketches: manysketch {sum(source.values())}, all-k zip {sum(every_k.values())}, k31 zip {sum(k31.values())}")
if every_k != source or k31 != want_k31:
    sys.exit("the re-written zips do not hold the same sketches as the manysketch zip")
PY
rm -rf "${d:?}"
echo "$(date -Is) REZIP_DONE"
