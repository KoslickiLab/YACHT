"""
Stage packaged YACHT databases for Zenodo: one folder per zip with the zip (hard link), metadata.json (Zenodo
deposition metadata), md5.txt and the data card, plus the scripts the database owner runs with their own token:
upload_to_zenodo.sh (creates draft depositions for the `yacht` community and uploads the files; publishes nothing),
update_zenodo_metadata.sh (gives the drafts the metadata of metadata.json after it changed) and
publish_zenodo_drafts.sh. They share zenodo_common.sh and zenodo_json.py.

Usage: stage_zenodo.py PACKAGE_DIR STAGING_DIR
"""
import json
import os
import re
import sys

MIN_YACHT_VERSION = "2.0.0"  # as in package_db.py

CREATORS = [{"name": "Koslicki, David", "affiliation": "Penn State University"}]
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


COMMON_SH = r"""# Shared by upload_to_zenodo.sh, update_zenodo_metadata.sh and publish_zenodo_drafts.sh (sourced after they cd here).
API=${ZENODO_API:-https://zenodo.org/api}
COMMUNITY=${ZENODO_COMMUNITY:-yacht}
: "${ZENODO_TOKEN:?set ZENODO_TOKEN to a Zenodo personal access token}"
JSON=(python3 zenodo_json.py)
WORK=$(mktemp -d .zenodo_tmp.XXXXXX)  # Zenodo's responses; in this folder, so that moving one into place is atomic
trap 'rm -rf "$WORK"' EXIT

# staged_dirs [prefix ...]: set DIRS to the given databases, or to every folder here (names may contain spaces)
staged_dirs() {
  local d
  DIRS=()
  if [ $# -gt 0 ]; then
    for d in "$@"; do DIRS+=("${d%/}"); done
  else
    for d in */; do
      if [ -d "$d" ]; then DIRS+=("${d%/}"); fi
    done
  fi
}

# zenodo METHOD URL OUT [curl option ...]: send a request to Zenodo, write the response body to OUT and print the
# HTTP status, or 000 if no complete response arrived. The token reaches curl in a header file, not on its command
# line, where other users of the machine could read it (ps).
zenodo() {
  local method=$1 url=$2 out=$3 code
  shift 3
  : > "$out"  # no stale answer if none arrives
  if code=$(curl -sS -X "$method" -o "$out" -w '%{http_code}' --connect-timeout 60 \
              -H @<(printf 'Authorization: Bearer %s\n' "$ZENODO_TOKEN") "$@" "$url"); then
    echo "$code"
  else
    echo 000
  fi
}

# status_text STATUS RESPONSE: "HTTP <status> (<Zenodo's message>)"
status_text() {
  local msg
  msg=$("${JSON[@]}" error "$2")
  if [ "$1" = 000 ]; then msg=${msg:-no complete answer}; fi
  echo "HTTP $1${msg:+ ($msg)}"
}

# check_auth STATUS RESPONSE: stop if Zenodo refused the token, since it would refuse every other request too
check_auth() {
  if [ "$1" = 401 ] || [ "$1" = 403 ]; then
    echo "Zenodo refused the request: $(status_text "$1" "$2"). Check the token and its scopes" \
         "(deposit:write, and deposit:actions to publish)." >&2
    exit 1
  fi
}

# keep RESPONSE DEST: move a checked response into place (a rename: DEST is either complete or unchanged)
keep() {
  mv -f "$1" "$2"
}
"""

JSON_PY = r'''"""
JSON checks for the Zenodo scripts of this folder (Zenodo's legacy deposition API). Each command prints its result,
or prints why it failed on stderr and exits with status 1.

  id DEPOSITION                        the id of the deposition that the file DEPOSITION describes
  get FILE KEY[.KEY...]                a value of the JSON object in FILE (not empty)
  url BASE NAME                        BASE/NAME, with NAME URL-encoded
  state DEPOSITION ID                  "draft", or "published <DOI> <URL>"; DEPOSITION must describe deposition ID
  files DEPOSITION MD5_TXT             checks that the deposition holds the zip of MD5_TXT, with that md5, and README.md
  new METADATA COMMUNITY               the request body that creates a draft with METADATA (metadata.json)
  update METADATA DEPOSITION COMMUNITY the request body that gives the draft DEPOSITION the metadata METADATA (the
                                       fields of METADATA and the community; the publication date is kept)
  differs METADATA DEPOSITION COMMUNITY  the metadata fields of DEPOSITION that differ from METADATA (none: empty)
  error RESPONSE                       Zenodo's error message in RESPONSE, if any (does not fail)
"""
import json
import sys
import urllib.parse


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(1)


def load(path):
    try:
        with open(path) as f:
            text = f.read()
    except OSError as e:
        fail(f"{path}: {e.strerror}")
    if not text.strip():
        fail(f"{path}: empty")
    try:
        data = json.loads(text)
    except ValueError:
        fail(f"{path}: not JSON")
    if not isinstance(data, dict):
        fail(f"{path}: not a JSON object")
    return data


def deposition_id(path):
    dep = load(path)
    dep_id = dep.get("id")
    if type(dep_id) is not int or dep_id <= 0:
        fail(f"{path}: no deposition id")
    return dep, dep_id


def wanted_metadata(metadata_path, community):
    metadata = load(metadata_path).get("metadata")
    if not isinstance(metadata, dict) or not metadata:
        fail(f"{metadata_path}: no metadata")
    return {**metadata, "communities": [{"identifier": community}]}


def comparable(key, value):
    # Zenodo returns each creator with all its fields, unset ones as null
    if key == "creators" and isinstance(value, list):
        return [{k: v for k, v in c.items() if v not in (None, "")} if isinstance(c, dict) else c for c in value]
    return value


def state(path, wanted_id):
    dep, dep_id = deposition_id(path)
    if str(dep_id) != wanted_id:
        fail(f"{path}: describes deposition {dep_id}, not {wanted_id}")
    submitted, status = dep.get("submitted"), dep.get("state")
    if submitted is True and status in ("done", "inprogress"):  # inprogress: published, and being edited
        doi = dep.get("doi") or (dep.get("metadata") or {}).get("doi")
        if not isinstance(doi, str) or not doi.startswith("10."):
            fail(f"{path}: deposition {dep_id} is published but has no DOI")
        links = dep.get("links") or {}
        return f"published {doi} {links.get('record_html') or links.get('html') or dep.get('doi_url') or ''}".strip()
    if submitted is False and status == "unsubmitted":
        return "draft"
    fail(f"{path}: deposition {dep_id} is neither a draft nor published (submitted: {submitted}, state: {status})")


def files(path, md5_txt):
    dep, dep_id = deposition_id(path)
    try:
        with open(md5_txt) as f:
            md5, name = f.read().strip().split(None, 1)
    except (OSError, ValueError):
        fail(f"{md5_txt}: expected '<md5>  <zip name>'")
    have = {f.get("filename"): f.get("checksum") for f in dep.get("files") or [] if isinstance(f, dict)}
    if have.get(name) != md5:
        fail(f"draft {dep_id} does not hold {name} with md5 {md5} (on Zenodo: {have.get(name) or 'missing'})")
    if "README.md" not in have:
        fail(f"draft {dep_id} does not hold README.md")


def differs(metadata_path, path, community):
    wanted = wanted_metadata(metadata_path, community)
    have = load(path).get("metadata") or {}
    return " ".join(k for k in wanted if comparable(k, wanted[k]) != comparable(k, have.get(k)))


def error(path):
    try:
        with open(path, errors="replace") as f:
            text = f.read()
    except OSError:
        return ""
    try:
        data = json.loads(text)
    except ValueError:
        return " ".join(text.split())[:300]
    if not isinstance(data, dict):
        return " ".join(text.split())[:300]
    parts = [str(data.get("message") or "")]
    for e in data.get("errors") or []:
        if isinstance(e, dict):
            parts.append(f"{e.get('field', '')}: {' '.join(map(str, e.get('messages') or [e.get('message', '')]))}")
    return "; ".join(p for p in parts if p)


def main():
    command, args = sys.argv[1], sys.argv[2:]
    if command == "id":
        print(deposition_id(args[0])[1])
    elif command == "get":
        value = load(args[0])
        for key in args[1].split("."):
            value = value.get(key) if isinstance(value, dict) else None
        if value in (None, ""):
            fail(f"{args[0]}: no {args[1]}")
        print(value)
    elif command == "url":
        print(f"{args[0]}/{urllib.parse.quote(args[1])}")
    elif command == "state":
        print(state(*args))
    elif command == "files":
        files(*args)
    elif command == "new":
        print(json.dumps({"metadata": wanted_metadata(*args)}))
    elif command == "update":
        metadata_path, path, community = args
        metadata = wanted_metadata(metadata_path, community)
        have = load(path).get("metadata") or {}
        if "publication_date" not in metadata and have.get("publication_date"):
            metadata["publication_date"] = have["publication_date"]  # Zenodo would set it to the day of the update
        print(json.dumps({"metadata": metadata}))
    elif command == "differs":
        print(differs(*args))
    elif command == "error":
        print(error(args[0]))
    else:
        fail(f"unknown command {command}")


if __name__ == "__main__":
    main()
'''

UPLOAD_SH = r"""#!/bin/bash
# Create one DRAFT Zenodo deposition per staged database and upload its zip and data card (README.md). The metadata
# (metadata.json) asks for inclusion in the `yacht` community: when a draft is published, the community owners are
# notified and accept it (on the community's "Requests" page, if it is not included automatically).
# Nothing is published here: review the drafts on https://zenodo.org/me/uploads, then publish them there or with
# publish_zenodo_drafts.sh. Re-running skips the databases whose upload is complete (zenodo_uploaded) and resumes
# the others in their existing draft (zenodo_deposition.json). The script stops at the first error.
#
# Usage: ZENODO_TOKEN=<token> bash upload_to_zenodo.sh [prefix ...]
#   token: https://zenodo.org/account/settings/applications/tokens/new/ with the scope deposit:write
#   (deposit:actions too, to publish with publish_zenodo_drafts.sh)
#   ZENODO_API=https://sandbox.zenodo.org/api (with a sandbox token) tries everything on the Zenodo sandbox first;
#   ZENODO_COMMUNITY=<id> uses another community (e.g. one created on the sandbox).
set -euo pipefail
cd "$(dirname "$0")"
source ./zenodo_common.sh
staged_dirs "$@"
for dir in "${DIRS[@]}"; do
  if [ -f "$dir/zenodo_uploaded" ]; then echo "== $dir: already uploaded, skipped"; continue; fi
  echo "== $dir"
  read -r md5 zip_name < "$dir/md5.txt"
  if [ ! -f "$dir/zenodo_deposition.json" ]; then  # otherwise: the draft of an interrupted run
    "${JSON[@]}" new "$dir/metadata.json" "$COMMUNITY" > "$WORK/body.json"
    code=$(zenodo POST "$API/deposit/depositions" "$WORK/deposition.json" -H "Content-Type: application/json" \
             --data-binary @"$WORK/body.json")
    check_auth "$code" "$WORK/deposition.json"
    if [[ $code != 2?? ]]; then
      echo "creating the draft failed: $(status_text "$code" "$WORK/deposition.json")"; exit 1
    fi
    "${JSON[@]}" id "$WORK/deposition.json" > /dev/null
    keep "$WORK/deposition.json" "$dir/zenodo_deposition.json"
  fi
  id=$("${JSON[@]}" id "$dir/zenodo_deposition.json")
  bucket=$("${JSON[@]}" get "$dir/zenodo_deposition.json" links.bucket)
  html=$("${JSON[@]}" get "$dir/zenodo_deposition.json" links.html)
  for f in "$dir/$zip_name" "$dir/README.md"; do
    url=$("${JSON[@]}" url "$bucket" "$(basename "$f")")
    code=$(zenodo PUT "$url" "$WORK/upload.json" --upload-file "$f")
    check_auth "$code" "$WORK/upload.json"
    if [[ $code != 2?? ]]; then echo "uploading $f failed: $(status_text "$code" "$WORK/upload.json")"; exit 1; fi
    if [ "$f" = "$dir/$zip_name" ]; then
      got=$("${JSON[@]}" get "$WORK/upload.json" checksum)
      [ "${got#md5:}" = "$md5" ] || { echo "md5 mismatch for $f: Zenodo has ${got#md5:}, expected $md5"; exit 1; }
    fi
  done
  touch "$dir/zenodo_uploaded"
  echo "draft $id $html (md5 $md5 verified)"
done
"""

UPDATE_SH = r"""#!/bin/bash
# Give the drafts on Zenodo the metadata of metadata.json, e.g. after the creators were changed in it. The fields of
# metadata.json and the community request are replaced; the publication date is kept. Zenodo's answer is compared
# with metadata.json, and it replaces zenodo_deposition.json. Drafts that already have this metadata are not changed,
# and published depositions are not changed (edit those on the website). Re-running is safe. Exit status 1 if a
# draft was not updated (they are listed at the end).
#
# Usage: ZENODO_TOKEN=<token> bash update_zenodo_metadata.sh [prefix ...]
#   token scope: deposit:write; ZENODO_API and ZENODO_COMMUNITY as for upload_to_zenodo.sh
set -euo pipefail
cd "$(dirname "$0")"
source ./zenodo_common.sh
staged_dirs "$@"

# update_one DIR (errors are handled here: set -e does not apply in a function called with ||)
update_one() {
  local dir=$1 id code state diff left
  if ! id=$("${JSON[@]}" id "$dir/zenodo_deposition.json"); then
    echo "== $dir: zenodo_deposition.json is missing or invalid, so the draft is unknown"; return 1
  fi
  code=$(zenodo GET "$API/deposit/depositions/$id" "$WORK/draft.json" --max-time 120)
  check_auth "$code" "$WORK/draft.json"
  if [ "$code" != 200 ]; then
    echo "== $dir: reading deposition $id failed: $(status_text "$code" "$WORK/draft.json")"; return 1
  fi
  if ! state=$("${JSON[@]}" state "$WORK/draft.json" "$id"); then echo "== $dir: not updated"; return 1; fi
  if [ "$state" != draft ]; then
    echo "== $dir: already published (${state#published }), not changed; edit it on the website"; return 1
  fi
  if ! diff=$("${JSON[@]}" differs "$dir/metadata.json" "$WORK/draft.json" "$COMMUNITY"); then
    echo "== $dir: not updated"; return 1
  fi
  if [ -z "$diff" ]; then
    keep "$WORK/draft.json" "$dir/zenodo_deposition.json" || return 1
    echo "== $dir: draft $id already has the metadata of metadata.json"; CURRENT=$((CURRENT + 1)); return 0
  fi
  if ! "${JSON[@]}" update "$dir/metadata.json" "$WORK/draft.json" "$COMMUNITY" > "$WORK/body.json"; then
    echo "== $dir: not updated"; return 1
  fi
  code=$(zenodo PUT "$API/deposit/depositions/$id" "$WORK/updated.json" -H "Content-Type: application/json" \
           --data-binary @"$WORK/body.json" --max-time 120)
  check_auth "$code" "$WORK/updated.json"
  if [ "$code" != 200 ]; then
    echo "== $dir: updating draft $id failed: $(status_text "$code" "$WORK/updated.json"); rerun"
    return 1
  fi
  if ! state=$("${JSON[@]}" state "$WORK/updated.json" "$id") || [ "$state" != draft ] ||
     ! left=$("${JSON[@]}" differs "$dir/metadata.json" "$WORK/updated.json" "$COMMUNITY") || [ -n "$left" ]; then
    echo "== $dir: Zenodo's answer to the update of draft $id is not this metadata (${left:-see above}); rerun"
    return 1
  fi
  keep "$WORK/updated.json" "$dir/zenodo_deposition.json" || return 1
  echo "== $dir: draft $id updated: $diff"; UPDATED=$((UPDATED + 1))
}

UPDATED=0; CURRENT=0; FAILED=()
for dir in "${DIRS[@]}"; do
  update_one "$dir" || FAILED+=("$dir")
done
echo "updated: $UPDATED; already up to date: $CURRENT; not updated: ${#FAILED[@]}"
if [ ${#FAILED[@]} -gt 0 ]; then printf '  %s\n' "${FAILED[@]}"; exit 1; fi
"""

PUBLISH_SH = r"""#!/bin/bash
# Publish the drafts made by upload_to_zenodo.sh. Publishing mints the DOIs and cannot be undone: review the drafts
# on https://zenodo.org/me/uploads first.
#
# For each database the deposition is read from Zenodo first. One that is already published (by an earlier run whose
# answer was lost, or on the website) is recorded, not published again. A draft is published only if its metadata
# on Zenodo is that of metadata.json (after changing metadata.json, run update_zenodo_metadata.sh) and it holds the
# zip (with the md5 of md5.txt) and README.md. zenodo_published.json is written only when Zenodo's answer describes
# this deposition as published, with its DOI; if the request fails or the answer is not that, the deposition is read
# again to find out whether it was published. Re-running is safe. Exit status 1 if a database is not published
# (they are listed at the end).
#
# Usage: ZENODO_TOKEN=<token> bash publish_zenodo_drafts.sh [prefix ...]
#   token scopes: deposit:write and deposit:actions; ZENODO_API and ZENODO_COMMUNITY as for upload_to_zenodo.sh
set -euo pipefail
cd "$(dirname "$0")"
source ./zenodo_common.sh
staged_dirs "$@"

# read_state ID: read deposition ID from Zenodo into $WORK/state.json and set STATE to "draft" or
# "published <DOI> <URL>" (up to 3 tries while Zenodo cannot be reached or answers 429/5xx)
read_state() {
  local code try
  STATE=
  for try in 1 2 3; do
    code=$(zenodo GET "$API/deposit/depositions/$1" "$WORK/state.json" --max-time 120)
    check_auth "$code" "$WORK/state.json"
    if [ "$code" = 200 ]; then
      STATE=$("${JSON[@]}" state "$WORK/state.json" "$1") && return 0
      return 1
    fi
    echo "   reading deposition $1: $(status_text "$code" "$WORK/state.json")"
    case $code in 000|429|5??) [ "$try" = 3 ] || sleep 10 ;; *) return 1 ;; esac
  done
  return 1
}

# record DIR RESPONSE: keep Zenodo's description of the published deposition as DIR/zenodo_published.json
record() {
  keep "$2" "$1/zenodo_published.json" ||
    { echo "== $1: published, but zenodo_published.json could not be written; rerun to record it"; return 1; }
}

# publish_one DIR (errors are handled here: set -e does not apply in a function called with ||)
publish_one() {
  local dir=$1 id code info diff what
  if [ ! -f "$dir/zenodo_uploaded" ]; then
    echo "== $dir: the upload is not complete (upload_to_zenodo.sh); not published"; return 1
  fi
  if ! id=$("${JSON[@]}" id "$dir/zenodo_deposition.json"); then
    echo "== $dir: zenodo_deposition.json is missing or invalid, so the draft is unknown; not published"; return 1
  fi
  if [ -e "$dir/zenodo_published.json" ]; then
    if info=$("${JSON[@]}" state "$dir/zenodo_published.json" "$id") && [ "${info%% *}" = published ]; then
      echo "== $dir: already published: ${info#published }"; ALREADY=$((ALREADY + 1)); return 0
    fi
    echo "   $dir/zenodo_published.json does not describe a published deposition; asking Zenodo"
  fi
  if ! read_state "$id"; then echo "== $dir: could not read deposition $id from Zenodo; not published"; return 1; fi
  if [ "${STATE%% *}" = published ]; then
    record "$dir" "$WORK/state.json" || return 1
    echo "== $dir: already published on Zenodo: ${STATE#published }"; ALREADY=$((ALREADY + 1)); return 0
  fi
  if ! diff=$("${JSON[@]}" differs "$dir/metadata.json" "$WORK/state.json" "$COMMUNITY"); then
    echo "== $dir: not published"; return 1
  fi
  if [ -n "$diff" ]; then
    echo "== $dir: draft $id on Zenodo differs from metadata.json in: $diff (run update_zenodo_metadata.sh);" \
         "not published"
    return 1
  fi
  if ! "${JSON[@]}" files "$WORK/state.json" "$dir/md5.txt"; then echo "== $dir: not published"; return 1; fi
  code=$(zenodo POST "$API/deposit/depositions/$id/actions/publish" "$WORK/publish.json" --max-time 600)
  check_auth "$code" "$WORK/publish.json"
  if [[ $code == 2?? ]]; then
    if info=$("${JSON[@]}" state "$WORK/publish.json" "$id" 2>&1) && [ "${info%% *}" = published ]; then
      record "$dir" "$WORK/publish.json" || return 1
      echo "== $dir: published ${info#published }"; PUBLISHED=$((PUBLISHED + 1)); return 0
    fi
    what="HTTP $code, but the answer does not describe it as published (${info#*: })"
  else
    what=$(status_text "$code" "$WORK/publish.json")
  fi
  echo "   publishing draft $id: $what; asking Zenodo whether it was published"
  if read_state "$id" && [ "${STATE%% *}" = published ]; then
    record "$dir" "$WORK/state.json" || return 1
    echo "== $dir: published ${STATE#published }"; PUBLISHED=$((PUBLISHED + 1)); return 0
  fi
  echo "== $dir: NOT published (deposition $id: ${STATE:-state unknown}; ${API%/api}/deposit/$id); rerun to retry," \
       "it records a publication it missed instead of publishing again"
  return 1
}

PUBLISHED=0; ALREADY=0; FAILED=()
for dir in "${DIRS[@]}"; do
  publish_one "$dir" || FAILED+=("$dir")
done
echo "published: $PUBLISHED; already published: $ALREADY; not published: ${#FAILED[@]}"
if [ ${#FAILED[@]} -gt 0 ]; then printf '  %s\n' "${FAILED[@]}"; exit 1; fi
"""

SCRIPTS = {"zenodo_common.sh": COMMON_SH, "zenodo_json.py": JSON_PY, "upload_to_zenodo.sh": UPLOAD_SH,
           "update_zenodo_metadata.sh": UPDATE_SH, "publish_zenodo_drafts.sh": PUBLISH_SH}


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
                "creators": CREATORS,
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
    for script, text in SCRIPTS.items():
        with open(os.path.join(staging, script), "w") as f:
            f.write(text)
    with open(os.path.join(staging, "INDEX.tsv"), "w") as f:
        f.write("prefix\tzip\tsize_bytes\tmd5\tlicense\n")
        f.writelines("\t".join(map(str, row)) + "\n" for row in staged)
    print(f"staged {len(staged)} databases in {staging}")


if __name__ == "__main__":
    main()
