#!/bin/bash
# Train one YACHT database and record time/memory.
# Usage: train_db.sh REF_ZIP PREFIX ANI OUTDIR [THREADS]
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
ref=$1; prefix=$2; ani=$3; outdir=$4; threads=${5:-64}
mkdir -p "$outdir"
echo "$(date -Is) TRAIN_START $prefix ani=$ani threads=$threads ref=$ref"
/usr/bin/time -v -o "$outdir/${prefix}_train_time.txt" env PATH="$YACHT_ENV/bin:$PATH" "$YACHT_ENV/bin/yacht" train \
  --ref_file "$ref" --ksize 31 --num_threads "$threads" --ani_thresh "$ani" --prefix "$prefix" --outdir "$outdir" --force \
  > "$outdir/${prefix}_train.log" 2>&1
rc=$?
echo "$(date -Is) TRAIN_END $prefix rc=$rc $(grep -E 'Elapsed|Maximum resident' "$outdir/${prefix}_train_time.txt" | awk -F': ' '{printf "%s=%s ", $1, $2}')"
