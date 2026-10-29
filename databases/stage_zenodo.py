"""
Stage packaged YACHT databases for Zenodo: one folder per zip with the zip (hard link), metadata.json (Zenodo
deposition metadata) and the data card, plus upload_to_zenodo.sh (run by the database owner with their own token;
it creates draft depositions in the `yacht` community and does not publish them).

Usage: stage_zenodo.py PACKAGE_DIR STAGING_DIR
"""
import json
import os
import re
import sys

LICENSES = {  # Zenodo license ids
    "gtdb": "cc-by-sa-4.0", "globdb": "cc-by-sa-4.0", "yacht-union": "cc-by-sa-4.0",
    "ictv": "cc-by-4.0", "genbank": "cc-zero", "ncbi-viruses": "cc-zero",
}
RELATED = {
    "gtdb": [("10.1093/nar/gkab776", "isDerivedFrom")],
    "globdb": [("10.1093/bioadv/vbaf280", "isDerivedFrom"), ("10.1093/nar/gkab776", "isDerivedFrom")],
    "yacht-union-globdb": [("10.1093/bioadv/vbaf280", "isDerivedFrom"), ("10.1093/nar/gkab776", "isDerivedFrom")],
    "yacht-union-insdc": [("10.1093/nar/gkab776", "isDerivedFrom")],
    "ictv": [("10.5281/zenodo.15042310", "isDerivedFrom")],
}
YACHT_PAPER = ("10.1093/bioinformatics/btae047", "isSupplementTo")


def family(prefix):
    for key in ("yacht-union-globdb", "yacht-union-insdc", "gtdb", "globdb", "ictv", "ncbi-viruses", "genbank"):
        if prefix.startswith(key):
            return key
    raise ValueError(prefix)


def main():
    package_dir, staging = sys.argv[1:3]
    os.makedirs(staging, exist_ok=True)
    staged = []
    for summary_path in sorted(f for f in os.listdir(package_dir) if f.endswith("_pretrained.summary.json")):
        summary = json.load(open(os.path.join(package_dir, summary_path)))
        zip_name = os.path.basename(summary["zip"])
        prefix = zip_name[: -len("_pretrained.zip")]
        fam = family(prefix)
        readme = open(os.path.join(package_dir, f"{prefix}_pretrained", "README.md")).read()
        title_line = readme.splitlines()[2] if len(readme.splitlines()) > 2 else prefix
        lic = next(v for k, v in LICENSES.items() if fam.startswith(k))
        related = [{"identifier": d, "relation": r, "scheme": "doi"} for d, r in RELATED.get(fam, []) + [YACHT_PAPER]]
        metadata = {
            "metadata": {
                "title": zip_name,
                "upload_type": "dataset",
                "description": (
                    f"<p>YACHT pretrained reference database <code>{prefix}</code>: {title_line}</p>"
                    f"<p>{summary['genomes']} genomes (from {summary['trained_on']} genome sketches; ANI threshold "
                    f"{re.sub(r'^.*_', '', prefix)}). Download with <code>yacht download pretrained_ref_db --name "
                    f"{prefix.rsplit('_', 1)[0]} --ani_thresh {prefix.rsplit('_', 1)[1]}</code> (YACHT 2.0 or later). "
                    f"Sources, licenses, identifiers and how to retrieve the genomes are in the README.md inside the zip.</p>"
                ),
                "creators": [{"name": "Koslicki, David"}],
                "license": lic,
                "keywords": ["YACHT", "metagenomics", "FracMinHash", "sourmash", fam],
                "communities": [{"identifier": "yacht"}],
                "related_identifiers": related,
                "version": "YACHT 2.0",
            }
        }
        out = os.path.join(staging, prefix)
        os.makedirs(out, exist_ok=True)
        target = os.path.join(out, zip_name)
        if not os.path.exists(target):
            os.link(summary["zip"], target)
        json.dump(metadata, open(os.path.join(out, "metadata.json"), "w"), indent=1)
        open(os.path.join(out, "README.md"), "w").write(readme)
        open(os.path.join(out, "md5.txt"), "w").write(f"{summary['md5']}  {zip_name}\n")
        staged.append((prefix, zip_name, summary["size"], summary["md5"], lic))
    with open(os.path.join(staging, "upload_to_zenodo.sh"), "w") as f:
        f.write('''#!/bin/bash
# Create one DRAFT Zenodo deposition per staged database (in the `yacht` community) and upload its zip.
# Nothing is published: review the drafts on https://zenodo.org/me/uploads and publish them there.
# Usage: ZENODO_TOKEN=<personal access token with deposit:write> bash upload_to_zenodo.sh [prefix ...]
set -euo pipefail
cd "$(dirname "$0")"
API=https://zenodo.org/api
for dir in ${@:-$(ls -d */ | tr -d /)}; do
  zip=$(ls $dir/*_pretrained.zip)
  echo "== $dir"
  dep=$(curl -sS -f -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" -X POST $API/deposit/depositions -d @$dir/metadata.json)
  bucket=$(echo "$dep" | python3 -c 'import json,sys; print(json.load(sys.stdin)["links"]["bucket"])')
  id=$(echo "$dep" | python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])')
  curl -sS -f -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file $zip "$bucket/$(basename $zip)" > /dev/null
  echo "draft deposition $id: https://zenodo.org/deposit/$id (md5 $(cut -d' ' -f1 $dir/md5.txt))"
done
''')
    with open(os.path.join(staging, "INDEX.tsv"), "w") as f:
        f.write("prefix\tzip\tsize_bytes\tmd5\tlicense\n")
        for row in staged:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"staged {len(staged)} databases in {staging}")


if __name__ == "__main__":
    main()
