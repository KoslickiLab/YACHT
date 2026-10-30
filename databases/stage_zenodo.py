"""
Stage packaged YACHT databases for Zenodo: one folder per zip with the zip (hard link), metadata.json (Zenodo
deposition metadata) and the data card, plus upload_to_zenodo.sh (run by the database owner with their own token;
it creates draft depositions for the `yacht` community and does not publish them) and publish_zenodo_drafts.sh.

Usage: stage_zenodo.py PACKAGE_DIR STAGING_DIR
"""
import json
import os
import re
import sys

MIN_YACHT_VERSION = "2.0.0"  # as in package_db.py

LICENSES = {  # Zenodo license ids
    "gtdb": "cc-by-sa-4.0", "globdb": "cc-by-sa-4.0", "yacht-union": "cc-by-sa-4.0",
    "ictv": "cc-by-4.0", "genbank": "cc-zero", "ncbi-viruses": "cc-zero",
}
RELATED = {
    "gtdb": [("10.1093/nar/gkab776", "isDerivedFrom")],
    "globdb": [("10.1093/bioadv/vbaf280", "isDerivedFrom"), ("10.1093/nar/gkab776", "isDerivedFrom")],
    "yacht-union-globdb": [("10.1093/bioadv/vbaf280", "isDerivedFrom"), ("10.1093/nar/gkab776", "isDerivedFrom"),
                           ("10.5281/zenodo.15042310", "isDerivedFrom")],
    "yacht-union-insdc": [("10.1093/nar/gkab776", "isDerivedFrom"), ("10.5281/zenodo.15042310", "isDerivedFrom")],
    "ictv": [("10.5281/zenodo.15042310", "isDerivedFrom")],
}
YACHT_PAPER = ("10.1093/bioinformatics/btae047", "isSupplementTo")


def family(prefix):
    for key in ("yacht-union-globdb", "yacht-union-insdc", "gtdb", "globdb", "ictv", "ncbi-viruses", "genbank"):
        if prefix.startswith(key):
            return key
    raise ValueError(prefix)


UPLOAD_SH = r"""#!/bin/bash
# Create one DRAFT Zenodo deposition per staged database and upload its zip and data card (README.md). The metadata
# (metadata.json) asks for inclusion in the `yacht` community: when a draft is published, the community owners are
# notified and accept it (on the community's "Requests" page, if it is not included automatically).
# Nothing is published here: review the drafts on https://zenodo.org/me/uploads, then publish them there or with
# publish_zenodo_drafts.sh. Re-running skips the databases whose upload is complete (zenodo_uploaded) and resumes
# the others in their existing draft (zenodo_deposition.json).
#
# Usage: ZENODO_TOKEN=<token> bash upload_to_zenodo.sh [prefix ...]
#   token: https://zenodo.org/account/settings/applications/tokens/new/ with the scope deposit:write
#   (deposit:actions too, to publish with publish_zenodo_drafts.sh)
#   ZENODO_API=https://sandbox.zenodo.org/api (with a sandbox token) tries everything on the Zenodo sandbox first;
#   ZENODO_COMMUNITY=<id> uses another community (e.g. one created on the sandbox).
set -euo pipefail
cd "$(dirname "$0")"
API=${ZENODO_API:-https://zenodo.org/api}
COMMUNITY=${ZENODO_COMMUNITY:-yacht}
: "${ZENODO_TOKEN:?set ZENODO_TOKEN to a Zenodo personal access token (scope deposit:write)}"
for dir in ${@:-$(ls -d */ | tr -d /)}; do
  if [ -f "$dir/zenodo_uploaded" ]; then echo "== $dir: already uploaded, skipped"; continue; fi
  zip=$(ls "$dir"/*_pretrained.zip)
  md5=$(cut -d' ' -f1 "$dir/md5.txt")
  echo "== $dir"
  if [ -f "$dir/zenodo_deposition.json" ]; then
    dep=$(cat "$dir/zenodo_deposition.json")  # draft of an interrupted run
  else
    meta=$(python3 -c 'import json,sys; m=json.load(open(sys.argv[1])); m["metadata"]["communities"]=[{"identifier": sys.argv[2]}]; print(json.dumps(m))' "$dir/metadata.json" "$COMMUNITY")
    dep=$(curl -sS -f -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" -X POST "$API/deposit/depositions" -d "$meta")
    echo "$dep" > "$dir/zenodo_deposition.json"
  fi
  bucket=$(echo "$dep" | python3 -c 'import json,sys; print(json.load(sys.stdin)["links"]["bucket"])')
  for f in "$zip" "$dir/README.md"; do
    res=$(curl -sS -f -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file "$f" "$bucket/$(basename "$f")")
    if [ "$f" = "$zip" ]; then
      got=$(echo "$res" | python3 -c 'import json,sys; print(json.load(sys.stdin)["checksum"].split(":")[-1])')
      [ "$got" = "$md5" ] || { echo "md5 mismatch for $zip: Zenodo has $got, expected $md5"; exit 1; }
    fi
  done
  touch "$dir/zenodo_uploaded"
  echo "draft $(echo "$dep" | python3 -c 'import json,sys; d=json.load(sys.stdin); print(d["id"], d["links"]["html"])') (md5 $md5 verified)"
done
"""

PUBLISH_SH = r"""#!/bin/bash
# Publish the drafts made by upload_to_zenodo.sh (zenodo_deposition.json in each folder). Publishing mints a DOI
# and cannot be undone: review the drafts first. Usage: ZENODO_TOKEN=<token with deposit:actions> bash
# publish_zenodo_drafts.sh [prefix ...]  (ZENODO_API as for upload_to_zenodo.sh)
set -euo pipefail
cd "$(dirname "$0")"
API=${ZENODO_API:-https://zenodo.org/api}
: "${ZENODO_TOKEN:?set ZENODO_TOKEN to a Zenodo personal access token (scopes deposit:write, deposit:actions)}"
for dir in ${@:-$(ls -d */ | tr -d /)}; do
  [ -f "$dir/zenodo_uploaded" ] || { echo "== $dir: no complete upload, skipped"; continue; }
  [ -f "$dir/zenodo_published.json" ] && { echo "== $dir: already published"; continue; }
  id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["id"])' "$dir/zenodo_deposition.json")
  curl -sS -f -H "Authorization: Bearer $ZENODO_TOKEN" -X POST "$API/deposit/depositions/$id/actions/publish" > "$dir/zenodo_published.json"
  echo "== $dir: published $(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d.get("doi"), d["links"].get("html"))' "$dir/zenodo_published.json")"
done
"""


def main():
    package_dir, staging = sys.argv[1:3]
    os.makedirs(staging, exist_ok=True)
    staged = []
    for summary_path in sorted(f for f in os.listdir(package_dir) if f.endswith("_pretrained.summary.json")):
        with open(os.path.join(package_dir, summary_path)) as f:
            summary = json.load(f)
        zip_name = os.path.basename(summary["zip"])
        prefix = zip_name[: -len("_pretrained.zip")]
        fam = family(prefix)
        with open(os.path.join(package_dir, f"{prefix}_pretrained", "README.md")) as f:
            readme = f.read()
        # the data card's title: the first paragraph after the heading that is not the version notice
        paragraphs = [p for p in readme.split("\n\n")[1:] if p.strip() and not p.startswith("**For YACHT")]
        title_line = " ".join(paragraphs[0].split()) if paragraphs else prefix
        lic = next(v for k, v in LICENSES.items() if fam.startswith(k))
        related = [{"identifier": d, "relation": r, "scheme": "doi"} for d, r in RELATED.get(fam, []) + [YACHT_PAPER]]
        metadata = {
            "metadata": {
                "title": f"{zip_name} (for YACHT >= {MIN_YACHT_VERSION})",
                "upload_type": "dataset",
                "description": (
                    f"<p><strong>For YACHT {MIN_YACHT_VERSION} or later; YACHT 1.x cannot use this database.</strong></p>"
                    f"<p>YACHT pretrained reference database <code>{prefix}</code>: {title_line}</p>"
                    f"<p>{summary['genomes']} genomes (from {summary['trained_on']} genome sketches; ANI threshold "
                    f"{re.sub(r'^.*_', '', prefix)}). Download with <code>yacht download pretrained_ref_db --name "
                    f"{prefix.rsplit('_', 1)[0]} --ani_thresh {prefix.rsplit('_', 1)[1]}</code>. "
                    f"Sources, licenses, identifiers, how coarse the database is and how to retrieve the genomes are in "
                    f"the data card README.md (in this record and inside the zip).</p>"
                ),
                "creators": [{"name": "Koslicki, David"}],
                "license": lic,
                "keywords": ["YACHT", f"YACHT {MIN_YACHT_VERSION}", "metagenomics", "FracMinHash", "sourmash", fam],
                "communities": [{"identifier": "yacht"}],
                "related_identifiers": related,
                "version": f"for YACHT >= {MIN_YACHT_VERSION}",
            }
        }
        out = os.path.join(staging, prefix)
        os.makedirs(out, exist_ok=True)
        target = os.path.join(out, zip_name)
        if not os.path.exists(target):
            os.link(summary["zip"], target)
        with open(os.path.join(out, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=1)
        with open(os.path.join(out, "README.md"), "w") as f:
            f.write(readme)
        with open(os.path.join(out, "md5.txt"), "w") as f:
            f.write(f"{summary['md5']}  {zip_name}\n")
        staged.append((prefix, zip_name, summary["size"], summary["md5"], lic))
    for script, text in (("upload_to_zenodo.sh", UPLOAD_SH), ("publish_zenodo_drafts.sh", PUBLISH_SH)):
        with open(os.path.join(staging, script), "w") as f:
            f.write(text)
    with open(os.path.join(staging, "INDEX.tsv"), "w") as f:
        f.write("prefix\tzip\tsize_bytes\tmd5\tlicense\n")
        f.writelines("\t".join(map(str, row)) + "\n" for row in staged)
    print(f"staged {len(staged)} databases in {staging}")


if __name__ == "__main__":
    main()
