#!/bin/bash
# Check a packaged database like a user would get it: unzip into a fresh folder, apply the steps of
# `yacht download pretrained_ref_db` after a download (config path rewrite, index check), then `yacht run` the SRA
# sample sketches (under /usr/bin/time; pinned to the core $CORE if set) and list the most-matched genomes that are
# found present.
# Usage: verify_db.sh ZIP SAMPLE_SKETCH...
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
zip=$1; shift
name=$(basename $zip .zip); W=$YACHT_DB_WORK/verify/$name
chmod -R u+rwx "${W:?}" 2>/dev/null; rm -rf "${W:?}"; mkdir -p $W; cp $zip $W/ && cd $W && unzip -q $(basename $zip) && rm $(basename $zip)
$YACHT_ENV/bin/python -I -c "
import sys
from yacht import download_pretrained_ref_db as d
d.update_config_file(sys.argv[1]); d.ensure_index(sys.argv[1])" $W/$name.zip > $W/postdownload.log 2>&1
cfg=$(ls $W/$name/*_config.json)
for s in "$@"; do
  acc=$(basename $s | cut -d. -f1); out=$W/run_$acc
  ${CORE:+taskset -c $CORE} /usr/bin/time -v -o $W/time_$acc.txt env PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/yacht run --json $cfg --sample_file $s --outdir $out > $W/run_$acc.log 2>&1
  rc=$?
  wall=$(grep Elapsed $W/time_$acc.txt | awk '{print $NF}'); rss=$(grep "Maximum resident" $W/time_$acc.txt | awk '{printf "%.1f", $NF/1e6}')
  $YACHT_ENV/bin/python -I - "$out/results/result_all.txt" "$acc" "$rc" "$wall" "$rss" <<'EOF'
import sys
import pandas as pd
path, acc, rc, wall, rss = sys.argv[1:6]
try:
    df = pd.read_csv(path, sep="\t")
except FileNotFoundError:
    print(f"{acc}: exit={rc} wall={wall} rss={rss}GB  (no result_all.txt)")
    sys.exit()
d = df[df["min_coverage"] == df["min_coverage"].min()]
present = d[d["in_sample_est"]]
top = present.sort_values("num_matches", ascending=False).head(3)
print(f"{acc}: exit={rc} wall={wall} rss={rss}GB overlapping={len(d)} present(min_cov={d['min_coverage'].min()})={len(present)}; top: "
      + " | ".join(f"{n[:60]} ({m})" for n, m in zip(top["organism_name"], top["num_matches"])))
EOF
done
