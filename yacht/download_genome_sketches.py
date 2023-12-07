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

BASE_URL = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/"

def generate_download_url(args):
    if args.database == "genbank":
        return f"{BASE_URL}{args.database}-{args.db_version}/{args.database}-{args.db_version}-{args.ncbi_organism}-k{args.k}.zip"
    elif args.database == "gtdb":
        suffix = "reps." if args.gtdb_type == "reps" else ""
        return f"{BASE_URL}{args.database}-{args.db_version}/{args.database}-{args.db_version}-{suffix}k{args.k}.zip"
    else:
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

def unzip_file(file_path, extract_to):
    subfolder_name = os.path.splitext(os.path.basename(file_path))[0]
    extract_path = os.path.join(extract_to, subfolder_name)

    if not os.path.exists(extract_path):
        logger.info(f"Creating subfolder for extraction: {extract_path}")
        os.makedirs(extract_path)

    logger.info(f"Starting to unzip {file_path} into {extract_path}")
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        logger.success(f"Successfully extracted {file_path} to {extract_path}")
    except zipfile.BadZipFile as e:
        logger.error(f"Failed to unzip {file_path}. Error: {e}")

def create_output_folder(outfolder):
    if not os.path.exists(outfolder):
        logger.info(f"Creating output folder: {outfolder}")
        os.makedirs(outfolder)

def main():
    parser = argparse.ArgumentParser(description="Download genome sketches for YACHT from the specified source.")
    parser.add_argument("--database", choices=['genbank', 'gtdb'], required=True)
    parser.add_argument("--db_version", required=True)
    parser.add_argument("--ncbi_organism", default="NULL")
    parser.add_argument("--gtdb_type", choices=[None, "reps", "full"], default=None)
    parser.add_argument("--k", choices=[21, 31, 51], type=int, required=True)
    parser.add_argument("--outfolder", help="Output folder for downloaded files.", default=".")

    args = parser.parse_args()

    download_url = generate_download_url(args)
    if not download_url:
        logger.error("Invalid URL generated from the given parameters.")
        return

    create_output_folder(args.outfolder)
    output_path = os.path.join(args.outfolder, os.path.basename(download_url))

    if download_file(download_url, output_path):
        unzip_file(output_path, args.outfolder)

if __name__ == "__main__":
    main()

