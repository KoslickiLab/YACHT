#!/usr/bin/env python
import requests
import argparse
from loguru import logger
import sys
import os
import json
import zipfile
import time
from tqdm import tqdm
from .utils import create_output_folder, check_download_args
# Import global variables
from .utils import ZENODO_COMMUNITY_URL

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
    parser.add_argument("--database", choices=["genbank", "gtdb"], required=True)
    parser.add_argument(
        "--db_version", choices=["genbank-2022.03", "rs214"], required=True
    )
    parser.add_argument(
        "--ncbi_organism",
        choices=["archaea", "bacteria", "fungi", "virus", "protozoa"],
        default=None,
    )
    parser.add_argument("--k", choices=[21, 31, 51], type=int, default=31)
    parser.add_argument(
        "--ani_thresh", type=float, choices=[0.80, 0.95, 0.995, 0.9995], required=True
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
    if args.database == "genbank":
        if args.db_version == "genbank-2022.03":
            return f"{args.db_version}-{args.ncbi_organism}-k{args.k}"
        else:
            logger.error(
                f"Invalid GenBank version: {args.db_version}. We currently only support genbank-2022.03."
            )
            return None
    else:
        if args.db_version == "rs214":
            return f"{args.database}-{args.db_version}-reps.k{args.k}"
        else:
            logger.error(
                f"Invalid GTDB version: {args.db_version}. We currently only support rs214."
            )
            return None


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
def unzip_file(file_path, extract_to):
    try:
        with zipfile.ZipFile(os.path.join(extract_to, file_path), 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        logger.success(f"Extracted {file_path} to {extract_to}")
    except zipfile.BadZipFile:
        logger.error(f"Failed to unzip {file_path}. It might not be a zip file.")


def main(args):
    ## Check if the input arguments are valid
    check_download_args(args, db_type="pretrained")

    ## Generate download URL
    file_prefix = generate_download_url(args)
    if not file_prefix:
        sys.exit(1)

    ## Fetch list of files from Zenodo
    zenodo_records = fetch_zenodo_records()
    if not zenodo_records:
        logger.error("No records fetched from Zenodo. Exiting.")
        sys.exit(1)
    current_pretrained_db_list = [record["title"] for record in zenodo_records]

    ## Create output folder if not exists
    create_output_folder(args.outfolder)

    ## Check if the specified db exists in Zenodo
    available_files = [
        filename for filename in current_pretrained_db_list if file_prefix in filename
    ]
    if len(available_files) == 0:
        logger.error(
            f"No pretrained database with prefix {file_prefix} found on Zenodo."
        )
        logger.info(f"Available pretrained databases: {current_pretrained_db_list}")
        sys.exit(1)
    else:
        available_files_with_ani = [
            filename for filename in available_files if f"{args.ani_thresh}" in filename
        ]
        available_ani_thresh = [float(x.split("_")[-2]) for x in available_files]
        if len(available_files_with_ani) == 0:
            logger.error(
                f"No pretrained database found for {file_prefix}_{args.ani_thresh}_pretrained.zip on Zenodo. "
                f"We currenlty only support {available_ani_thresh} for {file_prefix}."
            )
            sys.exit(1)
        else:
            file_name_to_search = available_files_with_ani[0]

    ## Check if the file already exists
    output_path = os.path.join(args.outfolder, file_name_to_search)
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        return

    ## Download the file
    file_url = next(
        (
            file.get("links", {}).get("self")
            for record in zenodo_records
            for file in record.get("files", [])
            if file_name_to_search in file.get("key", "")
        ),
        None,
    )

    if not file_url:
        logger.error(f"File '{file_name_to_search}' not found in Zenodo records.")
        sys.exit(1)
    
    if not download_file(file_url, output_path):
        logger.error(f"Failed to download '{file_name_to_search}' from Zenodo.")
        sys.exit(1)
    
    unzip_file(file_name_to_search, args.outfolder)
    update_config_file(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download pretrained models for YACHT from Zenodo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
