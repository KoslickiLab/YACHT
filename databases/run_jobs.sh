#!/bin/bash
# Run training jobs (lines "REF_ZIP PREFIX ANI") with at most P at a time. Usage: run_jobs.sh JOBS_FILE P THREADS
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
jobs=$1; P=$2; T=$3
export T
xargs -P $P -L 1 bash -c 'bash $SCRIPTS/train_db.sh "$0" "$1" "$2" $YACHT_DB_WORK/trained/${1}_pretrained $T' < $jobs
echo "$(date -Is) JOBS_DONE $jobs"
