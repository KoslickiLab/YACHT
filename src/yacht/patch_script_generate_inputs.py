#!/usr/bin/env python
import os
import sys
import argparse
import zipfile
from pathlib import Path
from loguru import logger
import json
import shutil
import pandas as pd
from tqdm import tqdm
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from utils import collect_signature_info, temp_generate_inputs

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

def add_arguments(parser):
    parser.add_argument(
        "--all_genome_name_path",
        type=str,
        help="Path to the file containing all the genome names.",
        required=True,
    )
    parser.add_argument(
        "--selected_genomes_file_path",
        type=str,
        help="Path to the file containing all the genome file path.",
        required=True,
    )
    parser.add_argument(
        "--sig_path_file",
        type=str,
        help="Path to the folder where the signature files are stored.",
        required=True,
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        help="Number of threads to use for parallelization.",
        required=False,
        default=16,
    )
    parser.add_argument(
        "--ksize",
        type=int,
        help="Size of kmers in sketch since Zipfiles can contain multiple k-sizes.",
        required=True,
    )
    parser.add_argument(
        "--ani_thresh",
        type=float,
        help="mutation cutoff for species equivalence. Organisms with this ANI "
        'or greater between them are considered "equivalent".',
        required=False,
        default=0.95,
    )
    parser.add_argument(
        "--prefix",
        help="Prefix name to identify this experiment.",
        required=False,
        default="yacht",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        help="Path to output directory.",
        required=False,
        default=os.getcwd(),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output directory if it exists.",
    )


def main(args):
    # get the arguments
    all_genome_name_path = str(Path(args.all_genome_name_path).absolute())
    selected_genomes_file_path = str(Path(args.selected_genomes_file_path).absolute())
    sig_path_file = str(Path(args.sig_path_file).absolute())
    num_threads = args.num_threads
    ksize = args.ksize
    ani_thresh = args.ani_thresh
    prefix = args.prefix
    outdir = str(Path(args.outdir).absolute())
    force = args.force

    # Create a temporary directory with time info as label
    logger.info("Creating a temporary directory")
    path_to_temp_dir = os.path.join(outdir, prefix + "_intermediate_files")
    if os.path.exists(path_to_temp_dir) and not force:
        raise ValueError(
            f"Temporary directory {path_to_temp_dir} already exists. Please remove it, use '--force', or given a new prefix name using parameter '--prefix'."
        )
    else:
        # remove the temporary directory if it exists
        if os.path.exists(path_to_temp_dir):
            logger.warning(
                f"Temporary directory {path_to_temp_dir} already exists."
            )
    os.makedirs(path_to_temp_dir, exist_ok=True)

    # Generate a folder `signatures` to store all the signature files
    signatures_folder = os.path.join(path_to_temp_dir, "signatures")
    os.makedirs(signatures_folder, exist_ok=True)
    
    # Copy the file from 'sig_path_file' to the new folder and rename it
    destination_file_path = os.path.join(path_to_temp_dir, "training_sig_files.txt")
    shutil.copy(sig_path_file, destination_file_path)
    
    # Create soft links for each *.sig file inside the 'signatures' folder
    logger.info("Create soft links for each *.sig file inside the 'signatures' folder")
    with open(destination_file_path, 'r') as file:
        for line in tqdm(file.readlines()):
            sig_path = line.strip()
            file_name = os.path.basename(sig_path)
            symlink_path = os.path.join(signatures_folder, file_name)
            if not os.path.exists(symlink_path):
                os.symlink(sig_path, symlink_path)

    # Extract signature information
    logger.info("Extracting signature information")
    sig_info_dict = collect_signature_info(num_threads, ksize, path_to_temp_dir)
    # check if all signatures have the same ksize and scaled
    logger.info("Checking if all signatures have the same scaled")
    scale_set = set([value[-2] for value in sig_info_dict.values()])
    if len(scale_set) != 1:
        raise ValueError(
            "Not all signatures have the same scaled. Please check your input."
        )
    scale = scale_set.pop()

    # Generate a dataframe with the selected genomes
    logger.info("Generating a dataframe with the selected genomes")
    manifest_df = temp_generate_inputs(selected_genomes_file_path, sig_info_dict, ksize, num_threads)

    # write out the manifest file
    logger.info("Writing out the manifest file")
    manifest_file_path = os.path.join(outdir, f"{prefix}_processed_manifest.tsv")
    manifest_df.to_csv(manifest_file_path, sep="\t", index=None)

    # save the config file
    logger.info("Saving the config file")
    json_file_path = os.path.join(outdir, f"{prefix}_config.json")
    json.dump(
        {
            "manifest_file_path": manifest_file_path,
            "intermediate_files_dir": path_to_temp_dir,
            "scale": scale,
            "ksize": ksize,
            "ani_thresh": ani_thresh,
        },
        open(json_file_path, "w"),
        indent=4,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
