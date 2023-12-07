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

def generate_file_name(args):
    if args.database == "genbank":
        return f"{args.database}-{args.db_version}-{args.ncbi_organism}-k{args.k}_{args.ani_thresh}_pretrained.zip"
    elif args.database == "gtdb":
        if args.gtdb_type == "full":
            return f"{args.database}-{args.db_version}-k{args.k}_{args.ani_thresh}_pretrained.zip"
        else:
            return f"{args.database}-{args.db_version}-{args.gtdb_type}.k{args.k}_{args.ani_thresh}_pretrained.zip"
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
    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        logger.success(f"Extracted {file_path} to {extract_to}")
    except zipfile.BadZipFile:
        logger.error(f"Failed to unzip {file_path}. It might not be a zip file.")

def create_output_folder(outfolder):
    if not os.path.exists(outfolder):
        logger.info(f"Creating output folder: {outfolder}")
        os.makedirs(outfolder)

def main():
    parser = argparse.ArgumentParser(description="Download pretrained models for YACHT from Zenodo.")
    parser.add_argument("--database", choices=['genbank', 'gtdb'], required=True)
    parser.add_argument("--db_version", required=True)
    parser.add_argument("--ncbi_organism", default=None)
    parser.add_argument("--gtdb_type", choices=[None, "reps", "full"], default=None)
    parser.add_argument("--k", choices=[21, 31, 51], type=int, required=True)
    parser.add_argument("--ani_thresh", choices=["0.80", "0.95", "0.995", "0.9995"], type=str, required=True)
    parser.add_argument("--outfolder", help="Output folder for downloaded files.", default=".")

    args = parser.parse_args()

    file_name_to_search = generate_file_name(args)
    if not file_name_to_search:
        logger.error("Invalid file name generated from the given parameters.")
        return

    zenodo_records = fetch_zenodo_records()
    if not zenodo_records:
        logger.error("No records fetched from Zenodo. Exiting.")
        return

    create_output_folder(args.outfolder)

    output_path = os.path.join(args.outfolder, file_name_to_search)
    if os.path.exists(output_path):
        logger.info(f"File {output_path} already exists. Skipping download.")
        unzip_file(output_path, args.outfolder)
        return

    file_url = next((file.get('links', {}).get('self') for record in zenodo_records for file in record.get('files', [])
                     if file_name_to_search in file.get('key', '')), None)

    if file_url and download_file(file_url, output_path):
        unzip_file(output_path, args.outfolder)
    else:
        logger.warning(f"File '{file_name_to_search}' not found in Zenodo records.")


if __name__ == "__main__":
    main()
