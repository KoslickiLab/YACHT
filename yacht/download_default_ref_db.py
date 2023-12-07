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
from .utils import BASE_URL

def generate_download_url(args):
    if args.database == "genbank":
        if args.db_version == "genbank-2022.03":
            if args.ncbi_organism == "virus":
                args.ncbi_organism = "viral"
            return f"{BASE_URL}{args.db_version}/{args.db_version}-{args.ncbi_organism}-k{args.k}.zip"
        else:
            logger.error(f"Invalid GenBank version: {args.db_version}. Now only support genbank-2022.03.")
            return None
    else:
        if args.db_version == "rs202":
            suffix = "-reps." if args.gtdb_type == "reps" else "."
            return f"{BASE_URL}{args.database}-{args.db_version}/{args.database}-{args.db_version}.genomic{suffix}k{args.k}.zip"
        elif args.db_version == "rs207":
            suffix = "-reps.dna." if args.gtdb_type == "reps" else "."
            return f"{BASE_URL}{args.database}-{args.db_version}/{args.database}-{args.db_version}.genomic{suffix}k{args.k}.zip"
        elif args.db_version == "rs214":
            suffix = "-reps." if args.gtdb_type == "reps" else "-"
            return f"{BASE_URL}{args.database}-{args.db_version}/{args.database}-{args.db_version}{suffix}k{args.k}.zip"
        else:
            logger.error(f"Invalid GTDB version: {args.db_version}. Now only support rs202, rs207, and rs214.")
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
    parser.add_argument("--db_version", choices=["genbank-2022.03", "rs202", "rs207", "rs214"], required=True)
    parser.add_argument("--ncbi_organism", choices=["archaea", "bacteria", "fungi", "virus", "protozoa"], default=None)
    parser.add_argument("--gtdb_type", choices=[None, "reps", "full"], default=None)
    parser.add_argument("--k", choices=[21, 31, 51], type=int, required=True)
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

    download_url = generate_download_url(args)
    if not download_url:
        os.exit(1)

    create_output_folder(args.outfolder)
    output_path = os.path.join(args.outfolder, os.path.basename(download_url))

    if download_file(download_url, output_path):
        unzip_file(output_path, args.outfolder)

if __name__ == "__main__":
    main()

