#!/usr/bin/env python
import requests
import argparse
from loguru import logger
import sys
import os
import zipfile

# Configure Loguru logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

# Import global variables
from .utils import ZENODO_COMMUNITY_URL

def fetch_zenodo_records():
    logger.info("Fetching list of files from Zenodo community 'yacht'")
    try:
        response = requests.get(ZENODO_COMMUNITY_URL)
        response.raise_for_status()
        return response.json().get('hits', {}).get('hits', [])
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching data from Zenodo: {e}")
        return []

def generate_download_url(args):
    if args.database == "genbank":
        if args.db_version == "genbank-2022.03":
            return f"{args.db_version}-{args.ncbi_organism}-k{args.k}"
        else:
            logger.error(f"Invalid GenBank version: {args.db_version}. Now only support genbank-2022.03.")
            return None
    else:
        if args.db_version == "rs214":
            return f"{args.database}-{args.db_version}-reps.k{args.k}"
        else:
            logger.error(f"Invalid GTDB version: {args.db_version}. Now only support rs214.")
            return None

def download_file(url, output_path):
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        return True
    try:
        logger.info(f"Starting download from {url}")
        response = requests.get(url)
        response.raise_for_status()
        with open(output_path, 'wb') as file:
            file.write(response.content)
        logger.success(f"Downloaded successfully and saved to {output_path}")
        return True
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download {url}: {e}")
        return False

def create_output_folder(outfolder):
    if not os.path.exists(outfolder):
        logger.info(f"Creating output folder: {outfolder}")
        os.makedirs(outfolder)

def main():
    parser = argparse.ArgumentParser(description="Download pretrained models for YACHT from Zenodo.")
    parser.add_argument("--database", choices=['genbank', 'gtdb'], required=True)
    parser.add_argument("--db_version", choices=["genbank-2022.03", "rs214"], required=True)
    parser.add_argument("--ncbi_organism", choices=["archaea", "bacteria", "fungi", "virus", "protozoa"], default=None)
    parser.add_argument("--k", choices=[21, 31, 51], type=int, required=True)
    parser.add_argument("--ani_thresh", type=float, choices=[0.80, 0.95, 0.995, 0.9995], required=True)
    parser.add_argument("--outfolder", help="Output folder for downloaded files.", default=".")

    args = parser.parse_args()

    ## Check if the input arguments are valid
    if args.database not in ["genbank", "gtdb"]:
        logger.error(f"Invalid database: {args.database}. Now only support genbank and gtdb.")
        os.exit(1)

    if args.k not in [21, 31, 51]:
        logger.error(f"Invalid k: {args.k}. Now only support 21, 31, and 51.")
        os.exit(1)

    if args.database == "genbank":
        if args.ncbi_organism is None:
            logger.warning("No NCBI organism specified using parameter --ncbi_organism. Use default: bacteria")
            args.ncbi_organism = "bacteria"
            
        if args.ncbi_organism not in ["archaea", "bacteria", "fungi", "virus", "protozoa"]:            
            logger.error(f"Invalid NCBI organism: {args.ncbi_organism}. Now only support archaea, bacteria, fungi, virus, and protozoa.")
            os.exit(1)
        
        if args.ncbi_organism == "virus":
            logger.error("We now have support for virus database.")
            os.exit(1)

    ## Generate download URL
    file_prefix = generate_download_url(args)
    if not file_prefix:
        os.exit(1)

    ## Fetch list of files from Zenodo
    zenodo_records = fetch_zenodo_records()
    if not zenodo_records:
        logger.error("No records fetched from Zenodo. Exiting.")
        os.exit(1)
    current_pretrained_db_list = [record['title'] for record in zenodo_records]

    ## Create output folder if not exists
    create_output_folder(args.outfolder)

    ## Check if the specified db exists in Zenodo
    available_files = [filename for filename in current_pretrained_db_list if file_prefix in filename]
    if len(available_files) == 0:
        logger.error(f"No pretrained database with prefix {file_prefix} found on Zenodo.")
        print(f"Available pretrained databases: {current_pretrained_db_list}")
        os.exit(1)
    else:
        available_files_with_ani = [filename for filename in available_files if f"{args.ani_thresh}" in filename]
        available_ani_thresh = [float(x.split('_')[-2]) for x in available_files]
        if len(available_files_with_ani) == 0:
            logger.error(f"No pretrained database found for {file_prefix}_{args.ani_thresh}_pretrained.zip on Zenodo. Now only support {available_ani_thresh} for {file_prefix}.")
            os.exit(1)
        else:
            file_name_to_search = available_files_with_ani[0]    

    ## Check if the file already exists
    output_path = os.path.join(args.outfolder, file_name_to_search)
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        return

    ## Download the file
    file_url = next((file.get('links', {}).get('self') for record in zenodo_records for file in record.get('files', [])
                     if file_name_to_search in file.get('key', '')), None)

    if file_url:
        download_file(file_url, output_path)
    else:
        logger.error(f"File {file_name_to_search} not found on Zenodo.")


if __name__ == "__main__":
    main()
