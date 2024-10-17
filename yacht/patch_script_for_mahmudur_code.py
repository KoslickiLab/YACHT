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
from utils import collect_signature_info, remove_corr_organisms_from_ref, extract_comparison_info

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
        "--comparion_path",
        type=str,
        help="Path to the folder where the individual comparison results are stored.",
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
    comparion_path = str(Path(args.comparion_path).absolute())
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
    scale_set = set([value[-1] for value in sig_info_dict.values()])
    if len(scale_set) != 1:
        raise ValueError(
            "Not all signatures have the same scaled. Please check your input."
        )
    scale = scale_set.pop()

    # Find the close related genomes with ANI > ani_thresh from the reference database
    logger.info(
        "Finding the closely related genomes with ANI > ani_thresh from the reference database"
    )
    comparison_result = extract_comparison_info(all_genome_name_path, comparion_path, ani_thresh, ksize, num_threads)

    # remove the close related organisms: any organisms with ANI > ani_thresh
    # pick only the one with largest number of unique kmers from all the closely related organisms
    logger.info("Removing the closely related organisms with ANI > ani_thresh")
    remove_corr_df, manifest_df = remove_corr_organisms_from_ref(
        sig_info_dict, comparison_result
    )

    # write out the manifest file
    logger.info("Writing out the manifest file")
    manifest_file_path = os.path.join(outdir, f"{prefix}_processed_manifest.tsv")
    manifest_df.to_csv(manifest_file_path, sep="\t", index=None)

    # write out a mapping dataframe from representative organism to the close related organisms
    logger.info(
        "Writing out a mapping dataframe from representative organism to the close related organisms"
    )
    if len(remove_corr_df) == 0:
        logger.warning("No close related organisms found.")
        remove_corr_df_indicator = ""
    else:
        remove_corr_df_path = os.path.join(
            outdir, f"{prefix}_removed_orgs_to_corr_orgas_mapping.tsv"
        )
        remove_corr_df.to_csv(remove_corr_df_path, sep="\t", index=None)
        remove_corr_df_indicator = remove_corr_df_path

    # save the config file
    logger.info("Saving the config file")
    json_file_path = os.path.join(outdir, f"{prefix}_config.json")
    json.dump(
        {
            "manifest_file_path": manifest_file_path,
            "remove_cor_df_path": remove_corr_df_indicator,
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
