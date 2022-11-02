#!/usr/bin/env bash
set -e
set -u
set -o pipefail

sourmash sketch dna -f -p k=31,scaled=1000,abund -o without_unknown_db.sig --singleton without_unknown_db.fasta

python ../../ref_matrix.py --ref_file without_unknown_db.sig  --ksize 31 --out_prefix default_EU
