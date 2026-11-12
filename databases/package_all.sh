#!/bin/bash
# Package every trained database (all four ANI thresholds) found in $YACHT_DB_WORK/trained. Usage: package_all.sh [name ...]
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
S=$YACHT_DB_WORK/src/sourmash; B=$YACHT_DB_WORK
declare -A KEY LIN TAX TAXLIN EXTRA
CONTAMINATED="$B/eukaryotes/fungi.ncbi-contaminated.csv $B/eukaryotes/protozoa.ncbi-contaminated.csv"
KEY[genbank-2025.04-fungi-k31]=genbank-2025.04-fungi;       LIN[genbank-2025.04-fungi-k31]="$S/genbank-2025.04/lineages.fungi.csv"
KEY[genbank-2025.04-protozoa-k31]=genbank-2025.04-protozoa; LIN[genbank-2025.04-protozoa-k31]="$S/genbank-2025.04/lineages.protozoa.csv"
KEY[ictv-vmr-msl40.k31-sc50]=ictv-vmr-msl40;                LIN[ictv-vmr-msl40.k31-sc50]="$S/2025-ictv-vmr-msl40/vmr_MSL40_v1.lineages.csv"
TAXLIN[ictv-vmr-msl40.k31-sc50]="$S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv"
KEY[ncbi-viruses-2025.01.k31-sc50]=ncbi-viruses-2025.01;    LIN[ncbi-viruses-2025.01.k31-sc50]="$S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv"
KEY[ncbi-viruses-refseq-2025.01.k31-sc50]=ncbi-viruses-refseq-2025.01; LIN[ncbi-viruses-refseq-2025.01.k31-sc50]="$S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv"
KEY[gtdb-rs226-reps.k31]=gtdb-rs226-reps;                   LIN[gtdb-rs226-reps.k31]="$S/gtdb-rs226/gtdb-rs226-reps.lineages.csv"
KEY[gtdb-rs232-reps.k31]=gtdb-rs232-reps;                   LIN[gtdb-rs232-reps.k31]="$B/gtdb-r232/gtdb-rs232-reps.lineages.csv"; TAX[gtdb-rs232-reps.k31]="$B/gtdb-r232/gtdb-rs232-reps.ncbi_taxids.tsv"
KEY[globdb-r232.k31]=globdb-r232;                           LIN[globdb-r232.k31]="$B/globdb-r232/globdb-r232.lineages.csv"; TAX[globdb-r232.k31]="$B/gtdb-r232/gtdb-rs232-reps.ncbi_taxids.tsv"
# unions: NCBI lineages for eukaryotes and RefSeq viruses, then the ICTV lineages (in the union's rank columns, written by
# assemble_unions.sh) for the VMR exemplars; NCBI taxids of the VMR exemplars from the NCBI viruses lineages
EUKVIR="$S/genbank-2025.04/lineages.fungi.csv $S/genbank-2025.04/lineages.protozoa.csv $S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv $B/unions/vmr_MSL40_v1.union-ranks.lineages.csv"
KEY[yacht-union-insdc-2026.k31]=yacht-union-insdc-2026;     LIN[yacht-union-insdc-2026.k31]="$B/gtdb-r232/gtdb-rs232-reps.lineages.csv $EUKVIR"; TAX[yacht-union-insdc-2026.k31]="$B/gtdb-r232/gtdb-rs232-reps.ncbi_taxids.tsv"
KEY[yacht-union-globdb-2026.k31]=yacht-union-globdb-2026;   LIN[yacht-union-globdb-2026.k31]="$B/globdb-r232/globdb-r232.lineages.csv $EUKVIR"; TAX[yacht-union-globdb-2026.k31]="$B/gtdb-r232/gtdb-rs232-reps.ncbi_taxids.tsv"
TAXLIN[yacht-union-insdc-2026.k31]="$S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv"; TAXLIN[yacht-union-globdb-2026.k31]="$S/ncbi-viruses-2025.01/ncbi-viruses.lineages.csv"
EXTRA[genbank-2025.04-fungi-k31]="$B/eukaryotes/fungi.ncbi-contaminated.csv"; EXTRA[genbank-2025.04-protozoa-k31]="$B/eukaryotes/protozoa.ncbi-contaminated.csv"
EXTRA[yacht-union-insdc-2026.k31]="$CONTAMINATED"; EXTRA[yacht-union-globdb-2026.k31]="$CONTAMINATED"
names=${@:-${!KEY[@]}}
$YACHT_ENV/bin/python -I $SCRIPTS/db_sources.py > /dev/null
for name in $names; do
  for ani in 0.80 0.95 0.995 0.9995; do
    prefix=${name}_$ani; d=$B/trained/${prefix}_pretrained
    [ -f $d/${prefix}_config.json ] || continue
    [ -f $B/package/${prefix}_pretrained.summary.json ] && continue
    tax=""; [ -n "${TAX[$name]:-}" ] && tax="--taxids ${TAX[$name]}"
    [ -n "${TAXLIN[$name]:-}" ] && tax="$tax --taxids-from-lineages ${TAXLIN[$name]}"
    [ -n "${EXTRA[$name]:-}" ] && tax="$tax --extra-files ${EXTRA[$name]}"
    echo "$(date -Is) package $prefix"
    PATH=$YACHT_ENV/bin:$PATH $YACHT_ENV/bin/python -I $SCRIPTS/package_db.py --db-dir $d --prefix $prefix --name $name --ani $ani \
      --sources $SCRIPTS/db_sources.json --key ${KEY[$name]} --lineages ${LIN[$name]} $tax --out-dir $B/package
  done
done
echo "$(date -Is) PACKAGE_DONE"
