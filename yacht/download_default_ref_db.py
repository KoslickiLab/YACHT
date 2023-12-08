#!/usr/bin/env python
import requests
import argparse
from loguru import logger
import sys
import os
from .utils import create_output_folder, check_download_args
import zipfile

# Configure Loguru logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

# Import global variables
from .utils import BASE_URL


def add_arguments(parser):
    parser.add_argument("--database", choices=['genbank', 'gtdb'], required=True)
    parser.add_argument("--db_version", choices=["genbank-2022.03", "rs202", "rs207", "rs214"], required=True)
    parser.add_argument("--ncbi_organism", choices=["archaea", "bacteria", "fungi", "virus", "protozoa"], default=None)
    parser.add_argument("--gtdb_type", choices=[None, "reps", "full"], default=None)
    parser.add_argument("--k", choices=[21, 31, 51], type=int, default=31)
    parser.add_argument("--outfolder", help="Output folder for downloaded files.", default=".")


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
        return True
    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download {url}: {e}")
        return False


def main(args):
    # Check if the input arguments are valid
    check_download_args(args)

    # Generate download URL
    download_url = generate_download_url(args)
    if not download_url:
        os.exit(1)

    # Create output folder if not exists
    create_output_folder(args.outfolder)
    output_path = os.path.join(args.outfolder, os.path.basename(download_url))

    # Download the file
    if download_file(download_url, output_path):
        logger.info(f"Downloaded successfully and saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download genome sketches for YACHT from the specified source.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args()
    main(args)

