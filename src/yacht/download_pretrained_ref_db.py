#!/usr/bin/env python
import requests
import argparse
from loguru import logger
import sys
import os
import json
import zipfile
import time
import re
from tqdm import tqdm
from .utils import create_output_folder, check_download_args
# Import global variables
from .utils import ZENODO_COMMUNITY_URL
from . import yacht_index

# Constants for retry logic
MAX_RETRIES = 3
DEFAULT_RETRY_WAIT = 20  # seconds

# Custom headers to avoid User-Agent based rate limiting
HEADERS = {'User-Agent': 'YACHT'}

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)




def add_arguments(parser):
    parser.add_argument(
        "--database",
        help="Database family, e.g. gtdb or genbank (see --list for all pretrained databases).",
    )
    parser.add_argument("--db_version", help="Database version, e.g. rs214 (gtdb) or genbank-2022.03 (genbank).")
    parser.add_argument(
        "--ncbi_organism",
        choices=["archaea", "bacteria", "fungi", "virus", "protozoa"],
        default=None,
    )
    parser.add_argument("--k", choices=[21, 31, 51], type=int, default=31)
    parser.add_argument(
        "--ani_thresh", type=float, help="ANI threshold the database was trained with (see --list)."
    )
    parser.add_argument(
        "--name",
        help="Name of a pretrained database as shown by --list (instead of --database, --db_version, --ncbi_organism and --k).",
    )
    parser.add_argument(
        "--list", action="store_true", help="List the pretrained databases available on Zenodo and exit."
    )
    parser.add_argument(
        "--outfolder", help="Output folder for downloaded files.", default="."
    )


def fetch_zenodo_records():
    logger.info("Fetching list of files from Zenodo community 'yacht'")
    all_records = []
    page = 1
    
    while True:
        url = f"{ZENODO_COMMUNITY_URL}&page={page}"
        
        # Retry logic for rate limiting
        for attempt in range(MAX_RETRIES):
            try:
                response = requests.get(url, headers=HEADERS)
                
                # Handle rate limiting (HTTP 429)
                if response.status_code == 429:
                    retry_after = int(response.headers.get('retry-after', DEFAULT_RETRY_WAIT))
                    logger.warning(f"Rate limited by Zenodo (429). Waiting {retry_after} seconds before retry (attempt {attempt + 1}/{MAX_RETRIES})...")
                    time.sleep(retry_after)
                    continue
                
                response.raise_for_status()
                data = response.json()
                hits = data.get("hits", {}).get("hits", [])
                if not hits:
                    logger.info(f"Fetched {len(all_records)} records from Zenodo")
                    return all_records
                all_records.extend(hits)
                # Check if we've fetched all records
                total = data.get("hits", {}).get("total", 0)
                if len(all_records) >= total:
                    logger.info(f"Fetched {len(all_records)} records from Zenodo")
                    return all_records
                page += 1
                break  # Success, move to next page
                
            except requests.exceptions.RequestException as e:
                if attempt < MAX_RETRIES - 1:
                    wait_time = DEFAULT_RETRY_WAIT * (attempt + 1)  # Exponential backoff
                    logger.warning(f"Request failed: {e}. Retrying in {wait_time} seconds (attempt {attempt + 1}/{MAX_RETRIES})...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"Error fetching data from Zenodo after {MAX_RETRIES} attempts: {e}")
                    # Return whatever records we collected so far
                    if all_records:
                        logger.warning(f"Returning {len(all_records)} records collected before failure")
                    return all_records
        else:
            # All retries exhausted for this page
            logger.error(f"Failed to fetch page {page} from Zenodo after {MAX_RETRIES} attempts")
            # Return whatever records we collected so far (page 1 might have succeeded)
            if all_records:
                logger.warning(f"Returning {len(all_records)} records collected before failure")
            return all_records
    
    return all_records


def generate_download_url(args):
    """
    The name (file name prefix) of the pretrained database selected by --database, --db_version, --ncbi_organism
    and --k: <db_version>-<ncbi_organism>-k<k> for genbank, gtdb-<db_version>-reps.k<k> for gtdb (the
    representative genomes), and <database>-<db_version>.k<k> for the other databases.
    """
    if args.database == "genbank":
        return f"{args.db_version}-{args.ncbi_organism}-k{args.k}"
    elif args.database == "gtdb":
        return f"{args.database}-{args.db_version}-reps.k{args.k}"
    else:
        return f"{args.database}-{args.db_version}.k{args.k}"


PRETRAINED_FILE_PATTERN = re.compile(r"^(?P<name>.+)_(?P<ani>[0-9]*\.?[0-9]+)_pretrained\.zip$")


def list_pretrained_databases(zenodo_records):
    """
    The pretrained databases in the Zenodo records: {name: {ani_thresh: (file name, download url, size, doi)}},
    from files named <name>_<ani_thresh>_pretrained.zip.
    """
    databases = {}
    for record in zenodo_records:
        for file in record.get("files", []):
            match = PRETRAINED_FILE_PATTERN.match(file.get("key", ""))
            if match:
                databases.setdefault(match["name"], {})[float(match["ani"])] = (
                    file["key"],
                    file.get("links", {}).get("self"),
                    file.get("size"),
                    record.get("doi"),
                )
    return databases


def print_pretrained_databases(databases):
    print("Pretrained databases available on Zenodo (https://zenodo.org/communities/yacht):")
    width = max((len(name) for name in databases), default=0)
    for name in sorted(databases):
        anis = ", ".join(f"{ani:g}" for ani in sorted(databases[name]))
        sizes = [entry[2] for entry in databases[name].values() if entry[2]]
        size = f"{min(sizes) / 1e9:.1f}-{max(sizes) / 1e9:.1f} GB" if sizes else ""
        print(f"  {name:<{width}}  ANI {anis}  {size}")
    print("Download one with: yacht download pretrained_ref_db --name <name> --ani_thresh <ANI> --outfolder <folder>")


def download_file(url, output_path):
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        return True
    
    for attempt in range(MAX_RETRIES):
        try:
            logger.info(f"Starting download from {url}")
            
            # Use streaming to avoid loading entire file into memory
            response = requests.get(url, stream=True, headers=HEADERS)
            
            # Handle rate limiting (HTTP 429)
            if response.status_code == 429:
                retry_after = int(response.headers.get('retry-after', DEFAULT_RETRY_WAIT))
                logger.warning(f"Rate limited by Zenodo (429). Waiting {retry_after} seconds before retry (attempt {attempt + 1}/{MAX_RETRIES})...")
                time.sleep(retry_after)
                continue
            
            response.raise_for_status()
            
            # Get file size for progress bar
            total_size = int(response.headers.get('content-length', 0))
            chunk_size = 8192  # 8KB chunks
            
            # Download with progress bar
            with open(output_path, "wb") as file:
                with tqdm(
                    total=total_size,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=os.path.basename(output_path),
                    disable=total_size == 0  # Disable if size unknown
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            file.write(chunk)
                            pbar.update(len(chunk))
            
            logger.success(f"Downloaded successfully and saved to {output_path}")
            return True
            
        except requests.exceptions.RequestException as e:
            if attempt < MAX_RETRIES - 1:
                wait_time = DEFAULT_RETRY_WAIT * (attempt + 1)  # Exponential backoff
                logger.warning(f"Download failed: {e}. Retrying in {wait_time} seconds (attempt {attempt + 1}/{MAX_RETRIES})...")
                time.sleep(wait_time)
            else:
                logger.error(f"Failed to download {url} after {MAX_RETRIES} attempts: {e}")
                return False
    
    return False

def update_config_file(file_path):
    try:
        absolute_path = os.path.abspath(file_path.replace(".zip",""))        
        config_file = [file for file in os.listdir(absolute_path) if "_config.json" in file][0]
        config_file = os.path.join(absolute_path,config_file)
        with open(config_file) as fp:
            config = json.loads(fp.read())
        for key in config:
            if isinstance(config[key], str) and config[key].startswith('/'):
                base_name = os.path.basename(config[key])
                config[key] = os.path.join(absolute_path, base_name)
        with open(config_file,'w') as fp:
            json.dump(config, fp, indent=4)
         
    except:
        logger.error(f"Could not find config file at {absolute_path}")
def ensure_index(file_path):
    """
    Make sure the extracted pretrained database has a valid YACHT index (needed by `yacht run`), creating it
    if the downloaded database does not include one. Never raises.
    :param file_path: path of the downloaded zip file (the database is extracted next to it)
    """
    folder = os.path.abspath(file_path.replace(".zip", ""))
    if os.path.isdir(folder):
        yacht_index.ensure_index_in_folder(folder)


def unzip_file(file_path, extract_to):
    try:
        with zipfile.ZipFile(os.path.join(extract_to, file_path), 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        logger.success(f"Extracted {file_path} to {extract_to}")
    except zipfile.BadZipFile:
        logger.error(f"Failed to unzip {file_path}. It might not be a zip file.")


def main(args):
    list_only = getattr(args, "list", False)
    name = getattr(args, "name", None)
    if not list_only:
        if args.ani_thresh is None:
            logger.error("Please give the ANI threshold of the database with --ani_thresh (see --list).")
            sys.exit(1)
        if name is None:
            if not args.database or not args.db_version:
                logger.error("Please select a database with --name, or with --database and --db_version (see --list).")
                sys.exit(1)
            ## Check if the input arguments are valid
            check_download_args(args, db_type="pretrained")

    ## Fetch list of files from Zenodo
    zenodo_records = fetch_zenodo_records()
    if not zenodo_records:
        logger.error("No records fetched from Zenodo. Exiting.")
        sys.exit(1)
    databases = list_pretrained_databases(zenodo_records)
    if list_only:
        print_pretrained_databases(databases)
        return

    ## Check if the specified db exists in Zenodo
    file_prefix = name if name is not None else generate_download_url(args)
    # e.g. --database ictv-vmr --db_version msl40 also selects ictv-vmr-msl40.k31-sc50
    candidates = [db for db in databases if db == file_prefix] or [
        db for db in databases if db.startswith(file_prefix + "-")
    ]
    if len(candidates) != 1:
        if candidates:
            logger.error(f"Several pretrained databases match {file_prefix}: {', '.join(sorted(candidates))}. Please choose one with --name.")
        else:
            logger.error(f"No pretrained database {file_prefix} found on Zenodo.")
            print_pretrained_databases(databases)
        sys.exit(1)
    file_prefix = candidates[0]
    available_ani_thresh = sorted(databases[file_prefix])
    matching_ani = [ani for ani in available_ani_thresh if abs(ani - args.ani_thresh) < 1e-9]
    if not matching_ani:
        logger.error(
            f"No pretrained database found for {file_prefix}_{args.ani_thresh}_pretrained.zip on Zenodo. "
            f"Available ANI thresholds for {file_prefix}: {', '.join(f'{ani:g}' for ani in available_ani_thresh)}."
        )
        sys.exit(1)
    file_name_to_search, file_url, _, _ = databases[file_prefix][matching_ani[0]]

    ## Create output folder if not exists
    create_output_folder(args.outfolder)

    ## Check if the file already exists
    output_path = os.path.join(args.outfolder, file_name_to_search)
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        # databases extracted by YACHT < 2.0 have no index yet
        ensure_index(output_path)
        return

    ## Download the file
    if not file_url:
        logger.error(f"File '{file_name_to_search}' not found in Zenodo records.")
        sys.exit(1)
    
    if not download_file(file_url, output_path):
        logger.error(f"Failed to download '{file_name_to_search}' from Zenodo.")
        sys.exit(1)
    
    unzip_file(file_name_to_search, args.outfolder)
    update_config_file(output_path)
    ensure_index(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download pretrained models for YACHT from Zenodo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
