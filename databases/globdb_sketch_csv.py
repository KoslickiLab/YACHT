"""
Build the `sourmash scripts manysketch` input CSV for the GlobDB r232 genomes without a GCA/GCF accession, their
lineages, and check that the GCA/GCF genomes of GlobDB r232 are exactly the GTDB R232 species representatives.

Usage: globdb_sketch_csv.py GLOBDB_RELEASE_DIR EXTRACTED_GENOMES_DIR GTDB_R232_LINEAGES_CSV OUT_CSV
  Writes OUT_CSV (name,genome_filename,protein_filename) with name = "<GlobDB ID> <GlobDB species>", and next to it
  <prefix>.lineages.csv (GlobDB taxonomy of the MAGs) and globdb-r232.lineages.csv (GTDB R232 part + MAGs).
"""
import csv
import gzip
import os
import re
import sys

release_dir, genomes_dir, gtdb_lineages, out_csv = sys.argv[1:5]
ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

taxonomy = {}
with gzip.open(os.path.join(release_dir, "globdb_r232_taxonomy.tsv.gz"), "rt") as f:
    for line in f:
        gid, lineage = line.rstrip("\n").split("\t")
        taxonomy[gid] = lineage.split(";")
accession_ids = {g for g in taxonomy if re.match(r"^GC[AF]_\d+$", g)}
mag_ids = set(taxonomy) - accession_ids
print(f"GlobDB r232 genomes: {len(taxonomy)}; with a GCA/GCF accession: {len(accession_ids)}; others: {len(mag_ids)}")

with open(gtdb_lineages, newline="") as f:
    gtdb_rows = list(csv.reader(f))
gtdb_unversioned = {r[0].rsplit(".", 1)[0] for r in gtdb_rows[1:]}
print(f"GTDB R232 representatives: {len(gtdb_unversioned)}; GlobDB accessions not in GTDB R232 reps: "
      f"{len(accession_ids - gtdb_unversioned)}; GTDB R232 reps not in GlobDB: {len(gtdb_unversioned - accession_ids)}")

files = {}
for root, _, names in os.walk(genomes_dir):
    for n in names:
        if n.endswith(".fa.gz"):
            files[n[: -len(".fa.gz")]] = os.path.join(root, n)
missing = sorted(mag_ids - set(files))
extra = sorted(set(files) - mag_ids)
print(f"MAG genome files: {len(files)}; MAGs without a file: {len(missing)} {missing[:5]}; files not in the taxonomy: {len(extra)} {extra[:5]}")

prefix = out_csv[: -len(".manysketch.csv")]
with open(out_csv, "w", newline="") as f_csv, open(prefix + ".lineages.csv", "w", newline="") as f_lin:
    w_csv, w_lin = csv.writer(f_csv), csv.writer(f_lin)
    w_csv.writerow(["name", "genome_filename", "protein_filename"])
    w_lin.writerow(["ident"] + ranks)
    for gid in sorted(mag_ids & set(files)):
        lineage = taxonomy[gid]
        assert len(lineage) == 7, (gid, lineage)
        species = re.sub(r"^s__", "", lineage[6])
        rest = " ".join(w for w in species.split() if w != gid)
        w_csv.writerow([f"{gid} {rest}" if rest else gid, files[gid], ""])
        w_lin.writerow([gid] + lineage)
with open(os.path.join(os.path.dirname(out_csv), "globdb-r232.lineages.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerows(gtdb_rows)
    with open(prefix + ".lineages.csv", newline="") as f_mag:
        w.writerows(list(csv.reader(f_mag))[1:])
print(f"wrote {out_csv}")
