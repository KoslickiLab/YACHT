#!/usr/bin/env python
import requests
import os
import sys
import argparse
from loguru import logger

# Configure Loguru logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")


GITHUB_API_URL = "https://api.github.com/repos/KoslickiLab/YACHT/contents/demo/{path}"
GITHUB_RAW_URL = "https://raw.githubusercontent.com/KoslickiLab/YACHT/main/demo/{path}"

def download_file(url, output_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'wb') as file:
            file.write(response.content)
        logger.info(f"Downloaded: {url}")
    else:
        logger.error(f"Failed to download {url}")

def fetch_file_list_from_github(folder_path=""):
    response = requests.get(GITHUB_API_URL.format(path=folder_path))
    if response.status_code == 200:
        file_list = response.json()
        return [file['path'].replace('demo/', '', 1) for file in file_list if file['type'] == 'file']
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

def main():
    parser = argparse.ArgumentParser(description="Download YACHT demo files.")
    parser.add_argument("--output", help="Output folder.", default="demo")
    args = parser.parse_args()

    logger.info(f"Starting download of YACHT demo files to {args.output}")
    download_demo_files(args.output)
    logger.info("Download completed.")

if __name__ == "__main__":
    main()

