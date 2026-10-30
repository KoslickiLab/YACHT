"""
Package a database trained with `yacht train` as a YACHT 2.0 pretrained database zip for Zenodo.

Usage: package_db.py --db-dir DIR --prefix PREFIX --name NAME --ani ANI --sources SOURCES_JSON --key KEY
                     --lineages CSV [--lineages CSV ...] [--taxids TSV ...] --out-dir OUT

The zip <PREFIX>_pretrained.zip holds the folder <PREFIX>_pretrained/ (the layout `yacht download
pretrained_ref_db` expects) with:
  <PREFIX>_config.json, <PREFIX>_processed_manifest.tsv, <PREFIX>_index.yidx   what `yacht run` needs
  <PREFIX>_lineages.csv          sourmash lineage CSV of the genomes in the database (ident = first word of the name)
  <PREFIX>_genome_to_taxid.tsv   genome_id (= organism_name in YACHT's output) and NCBI taxid, for `yacht convert`
  <PREFIX>_removed_genomes.tsv   genomes left out by `yacht train` and the genomes they were found equivalent to
  README.md                      data card
The decompressed reference signatures (<PREFIX>_intermediate_files) are not included: `yacht run` reads the index.
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import zipfile
from multiprocessing import Pool

import pandas as pd

csv.field_size_limit(sys.maxsize)


def ident(name):
    return name.split()[0] if name.split() else name


def read_lineages(paths):
    """{ident: {rank: value}} from sourmash lineage CSVs (ident column 'ident' or 'identifiers')."""
    lineages, ranks = {}, []
    for path in paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            key = "ident" if "ident" in reader.fieldnames else "identifiers"
            for r in reader.fieldnames:
                if r not in (key, "taxid", "taxpath", "name") and r not in ranks:
                    ranks.append(r)
            for row in reader:
                lineages[row[key]] = row
    return lineages, ranks


def read_taxids(paths, lineages, taxid_lineages=()):
    taxids = {}
    for path in taxid_lineages:  # only the taxid column of other lineage CSVs (same versioned accession)
        extra, _ = read_lineages([path])
        taxids.update({i: row["taxid"] for i, row in extra.items() if row.get("taxid")})
    taxids.update({i: row["taxid"] for i, row in lineages.items() if row.get("taxid")})
    for path in paths:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        taxids.update(dict(zip(df["ident"], df["ncbi_taxid"])))
    return {i: t for i, t in taxids.items() if t and t != "none"}


def _closest_kept(args):
    """removed genome index -> (highest containment in a kept genome, that kept genome, number of kept genomes it is
    equivalent to), from one comparison file (lines i,j,jaccard,containment of i in j,containment of j in i)"""
    path, kept_index = args
    out = {}
    with open(path) as f:
        for line in f:
            i, j, _, c = line.split(",", 4)[:4]
            i, j = int(i), int(j)
            if j in kept_index and i not in kept_index:
                c = float(c)
                best = out.get(i)
                if best is None:
                    out[i] = (c, j, 1)
                elif (c, -j) > (best[0], -best[1]):  # highest containment; ties: smallest genome index
                    out[i] = (c, j, best[2] + 1)
                else:
                    out[i] = (best[0], best[1], best[2] + 1)
    return out


def removed_genomes(intermediate_dir, kept_names, ksize, previous, processes=32):
    """Rows (removed genome id, name, closest genome in the database, containment in it, ANI estimate, number of genomes
    in the database it is equivalent to). When the comparison files of yacht train are gone (deleted after an earlier
    packaging), the rows of that packaging (the TSV `previous`) are used."""
    if not os.path.isdir(os.path.join(intermediate_dir, "comparison_files")) and os.path.exists(previous):
        with open(previous, newline="") as f:
            rows = [tuple(r) for r in csv.reader(f, delimiter="\t")][1:]
        return rows, len(kept_names) + len(rows)
    with open(os.path.join(intermediate_dir, "SOURMASH-MANIFEST.csv")) as f:
        next(f)
        # <md5>.sig.gz -> <md5>.sig and <md5>.sig.gz_<n> -> <md5>.sig_<n>, as yacht train decompresses them
        mf = {re.sub(r"\.gz(_\d+)?$", lambda m: m.group(1) or "", os.path.basename(r["internal_location"])): r
              for r in csv.DictReader(f)}
    with open(os.path.join(intermediate_dir, "training_sig_files.tsv")) as f:
        names = [mf[os.path.basename(line.strip().strip('"'))]["name"] for line in f if line.strip()]
    kept = set(kept_names)
    kept_index = frozenset(i for i, n in enumerate(names) if n in kept)
    files = glob.glob(os.path.join(intermediate_dir, "comparison_files", "*.txt"))
    closest = {}
    with Pool(processes) as pool:
        for part in pool.imap_unordered(_closest_kept, [(p, kept_index) for p in files]):
            for i, (c, j, n) in part.items():
                old = closest.get(i)
                if old is None:
                    closest[i] = (c, j, n)
                else:
                    best = (c, j) if (c, -j) > (old[0], -old[1]) else (old[0], old[1])
                    closest[i] = (best[0], best[1], old[2] + n)
    rows = []
    for i, name in enumerate(names):
        if i in kept_index:
            continue
        c, j, n = closest.get(i, (None, None, 0))
        rows.append((ident(name), name, ident(names[j]) if j is not None else "",
                     f"{c:.6g}" if c is not None else "", f"{c ** (1 / ksize):.6g}" if c is not None else "", n))
    return rows, len(names)


def coarseness(lineages, ranks, kept_ids, rows):
    """Markdown table: per rank, the taxa of all genomes yacht train was trained on, those that keep a genome in the
    database, and the share of the left-out genomes whose closest database genome is in another taxon at that rank or
    above (ranks without a value in either lineage are skipped); and the number of left-out genomes counted. NCBI's
    placeholders for missing ranks ("unclassified Caudoviricetes family") count as no value."""
    def known(lineage):
        return {k: v for k, v in lineage.items() if v and not v.startswith("unclassified ")}

    trained = [i for i in list(kept_ids) + [r[0] for r in rows] if i in lineages]
    kept = [i for i in kept_ids if i in lineages]
    paired = [(known(lineages[r[0]]), known(lineages[r[2]])) for r in rows if r[2] and r[0] in lineages and r[2] in lineages]
    lineages = {i: known(lineages[i]) for i in trained}
    differs = [False] * len(paired)
    lines = ["| rank | taxa | taxa with a genome in this database | left-out genomes represented by a genome of another taxon |",
             "|---|---:|---:|---:|"]
    for rank in ranks:
        for k, (a, b) in enumerate(paired):
            if a.get(rank) and b.get(rank) and a[rank] != b[rank]:
                differs[k] = True
        taxa = {lineages[i][rank] for i in trained if lineages[i].get(rank)}
        if not taxa:
            continue
        kept_taxa = {lineages[i][rank] for i in kept if lineages[i].get(rank)}
        share = f"{100 * sum(differs) / len(paired):.1f}%" if paired else "-"
        lines.append(f"| {rank} | {len(taxa):,} | {len(kept_taxa):,} ({100 * len(kept_taxa) / len(taxa):.1f}%) | {share} |")
    return "\n".join(lines), len(paired)


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--db-dir", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--ani", required=True)
    p.add_argument("--sources", required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--lineages", nargs="+", required=True)
    p.add_argument("--taxids", nargs="*", default=[])
    p.add_argument("--taxids-from-lineages", nargs="*", default=[])
    p.add_argument("--extra-files", nargs="*", default=[], help="files added to the zip as they are (e.g. lists of left-out genomes)")
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()

    with open(a.sources) as f:
        sources = json.load(f)
    card = sources[a.key]
    config_path = os.path.join(a.db_dir, f"{a.prefix}_config.json")
    with open(config_path) as f:
        config = json.load(f)
    manifest = pd.read_csv(config["manifest_file_path"], sep="\t", header=0)
    kept_names = manifest["organism_name"].astype(str).tolist()

    folder_name = f"{a.prefix}_pretrained"
    stage = os.path.join(a.out_dir, folder_name)
    os.makedirs(stage, exist_ok=True)

    # lineages and taxids of the genomes in the database
    lineages, ranks = read_lineages(a.lineages)
    taxids = read_taxids(a.taxids, lineages, a.taxids_from_lineages)
    lin_path = os.path.join(stage, f"{a.prefix}_lineages.csv")
    tax_path = os.path.join(stage, f"{a.prefix}_genome_to_taxid.tsv")
    no_lineage, no_taxid = [], []
    with open(lin_path, "w", newline="") as f_lin, open(tax_path, "w", newline="") as f_tax:
        w_lin = csv.writer(f_lin)
        w_tax = csv.writer(f_tax, delimiter="\t", lineterminator="\n")
        w_lin.writerow(["ident"] + ranks)
        w_tax.writerow(["genome_id", "taxid"])
        for name in kept_names:
            i = ident(name)
            if i in lineages:
                w_lin.writerow([i] + [lineages[i].get(r, "") for r in ranks])
            else:
                no_lineage.append(i)
            if i in taxids:
                w_tax.writerow([name, taxids[i]])
            else:
                no_taxid.append(i)

    # genomes left out by yacht train
    rem_path = os.path.join(stage, f"{a.prefix}_removed_genomes.tsv")
    rows, n_trained = removed_genomes(config["intermediate_files_dir"], kept_names, config["ksize"], rem_path)
    with open(rem_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["removed_genome", "removed_genome_name", "closest_genome_in_database", "containment_in_closest",
                    "ani_estimate", "n_equivalent_genomes_in_database"])
        w.writerows(rows)

    # config: the paths are rewritten to the extracted folder by `yacht download pretrained_ref_db`
    out_config = dict(config)
    for k in ("manifest_file_path", "intermediate_files_dir", "index_file_path"):
        out_config[k] = os.path.join("/", folder_name, os.path.basename(config[k]))
    with open(os.path.join(stage, os.path.basename(config_path)), "w") as f:
        json.dump(out_config, f, indent=4)

    yacht_version = subprocess.run([sys.executable, "-c", "import yacht; print(yacht.__version__)"],
                                   capture_output=True, text=True, check=True).stdout.strip()
    table, n_paired = coarseness(lineages, ranks, [ident(n) for n in kept_names], rows)
    warning = ""
    if float(a.ani) <= 0.80:
        warning = sources["ani_0.80_warning"].format(ksize=config["ksize"], containment=f"{float(a.ani) ** config['ksize']:.2g}")
    card_text = card["readme"].format(
        prefix=a.prefix, name=a.name, ani=a.ani, n_genomes=len(kept_names), n_trained=n_trained,
        n_removed=len(rows), ksize=config["ksize"], scaled=config["scale"], yacht_version=yacht_version,
        n_no_lineage=len(no_lineage), n_no_taxid=len(no_taxid), ani_warning=warning, coarseness=table,
        n_paired=n_paired,
    )
    with open(os.path.join(stage, "README.md"), "w") as f:
        f.write(card_text)

    zip_path = os.path.join(a.out_dir, f"{folder_name}.zip")
    files = [
        (config["manifest_file_path"], True), (config["index_file_path"], False),
        (os.path.join(stage, os.path.basename(config_path)), True), (lin_path, True), (tax_path, True),
        (rem_path, True), (os.path.join(stage, "README.md"), True),
    ] + [(path, True) for path in a.extra_files]
    with zipfile.ZipFile(zip_path, "w", allowZip64=True) as z:
        folder = zipfile.ZipInfo(folder_name + "/")
        folder.external_attr = (0o40755 << 16) | 0x10  # drwxr-xr-x (without it, `unzip` creates the folder with mode 000)
        z.writestr(folder, "")
        for path, compress in files:
            z.write(path, f"{folder_name}/{os.path.basename(path)}",
                    compress_type=zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED)
    summary = {
        "zip": zip_path, "md5": md5sum(zip_path), "size": os.path.getsize(zip_path), "genomes": len(kept_names),
        "trained_on": n_trained, "removed": len(rows), "without_lineage": len(no_lineage),
        "without_ncbi_taxid": len(no_taxid), "examples_without_lineage": no_lineage[:5],
    }
    with open(os.path.join(a.out_dir, f"{folder_name}.summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
