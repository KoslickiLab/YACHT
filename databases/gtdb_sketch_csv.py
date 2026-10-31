"""
Build the `sourmash scripts manysketch` input CSV for the GTDB R232 species representatives, plus their lineages
(sourmash lineage CSV) and NCBI taxids.

Usage: gtdb_sketch_csv.py GTDB_RELEASE_DIR EXTRACTED_GENOMES_DIR OUT_CSV
  Writes OUT_CSV (name,genome_filename,protein_filename) with name = "<versioned accession> <NCBI organism name>",
  and next to it <prefix>.lineages.csv (ident,superkingdom,phylum,class,order,family,genus,species) and
  <prefix>.ncbi_taxids.tsv (ident, ncbi_taxid, ncbi_species_taxid), ident = versioned accession.
"""
import csv
import os
import sys

import pandas as pd

release_dir, genomes_dir, out_csv = sys.argv[1:4]
cols = ["accession", "gtdb_representative", "gtdb_taxonomy", "ncbi_organism_name", "ncbi_taxid", "ncbi_species_taxid"]
meta = pd.concat(
    pd.read_csv(os.path.join(release_dir, f), sep="\t", usecols=cols, dtype=str, keep_default_na=False)
    for f in ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz")
)
reps = meta[meta["gtdb_representative"] == "t"].copy()
reps["ident"] = reps["accession"].str.replace(r"^(RS_|GB_)", "", regex=True)
assert reps["ident"].is_unique
by_ident = reps.set_index("ident")
print(f"species representatives in the metadata: {len(reps)}")

files = {}
for root, _, names in os.walk(genomes_dir):
    for n in names:
        if n.endswith("_genomic.fna.gz"):
            files[n[: -len("_genomic.fna.gz")]] = os.path.join(root, n)
missing = sorted(set(by_ident.index) - set(files))
extra = sorted(set(files) - set(by_ident.index))
print(f"genome files: {len(files)}; representatives without a file: {len(missing)} {missing[:5]}; files that are not representatives: {len(extra)} {extra[:5]}")

prefix = out_csv[: -len(".manysketch.csv")] if out_csv.endswith(".manysketch.csv") else os.path.splitext(out_csv)[0]
ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
with open(out_csv, "w", newline="") as f_csv, open(prefix + ".lineages.csv", "w", newline="") as f_lin, \
        open(prefix + ".ncbi_taxids.tsv", "w", newline="") as f_tax:
    w_csv, w_lin, w_tax = csv.writer(f_csv), csv.writer(f_lin), csv.writer(f_tax, delimiter="\t")
    w_csv.writerow(["name", "genome_filename", "protein_filename"])
    w_lin.writerow(["ident"] + ranks)
    w_tax.writerow(["ident", "ncbi_taxid", "ncbi_species_taxid"])
    for ident in sorted(set(files) & set(by_ident.index)):
        row = by_ident.loc[ident]
        organism = row["ncbi_organism_name"].strip()
        w_csv.writerow([f"{ident} {organism}" if organism else ident, files[ident], ""])
        lineage = row["gtdb_taxonomy"].split(";")
        assert len(lineage) == 7, (ident, lineage)
        w_lin.writerow([ident] + lineage)
        w_tax.writerow([ident, row["ncbi_taxid"], row["ncbi_species_taxid"]])
print(f"wrote {out_csv}")
