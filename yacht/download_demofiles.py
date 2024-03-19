#!/usr/bin/env python
import requests
import os
import sys
import argparse
from loguru import logger
# Import global variables
from .utils import GITHUB_API_URL, GITHUB_RAW_URL

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)




def add_arguments(parser):
    parser.add_argument("--outfolder", help="Output folder.", default="demo")


def download_file(url, output_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            file.write(response.content)
        logger.info(f"Downloaded: {url}")
    else:
        logger.error(f"Failed to download {url}")


def fetch_file_list_from_github(folder_path=""):
    response = requests.get(GITHUB_API_URL.format(path=folder_path))
    if response.status_code == 200:
        file_list = response.json()
        return [
            file["path"].replace("demo/", "", 1)
            for file in file_list
            if file["type"] == "file"
        ]
    else:
        logger.error(f"Failed to fetch file list from {folder_path}")
        return []


def download_demo_files(output_folder):
    folders_to_download = ["", "query_data", "ref_genomes"]
    for folder in folders_to_download:
        files_to_download = fetch_file_list_from_github(folder)
        for file_path in files_to_download:
            output_file_path = os.path.join(output_folder, file_path)
            os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
            download_file(GITHUB_RAW_URL.format(path=file_path), output_file_path)


def main(args):
    logger.info(f"Starting download of YACHT demo files to {args.outfolder}")
    download_demo_files(args.outfolder)
    logger.info("Download completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download YACHT demo files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
