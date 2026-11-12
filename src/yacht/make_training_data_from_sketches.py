#!/usr/bin/env python
import os
import sys
import argparse
import zipfile
from pathlib import Path
from loguru import logger
import json
import shutil
import glob
from . import utils
from .utils import __version__ as YACHT_VERSION
from .yacht_index import INDEX_SUFFIX

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)


def _positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1 (got {value})")
    return number


def _ani_threshold(value):
    number = float(value)
    if not (0 < number <= 1):
        raise argparse.ArgumentTypeError(f"must be greater than 0 and at most 1 (got {value})")
    return number


def add_arguments(parser):
    parser.add_argument(
        "--ref_file",
        help="Location of the Sourmash signature database file. "
        "This is expected to be in Zipfile format (eg. *.zip) "
        'that contains a manifest "SOURMASH-MANIFEST.csv" and a folder "signatures" '
        "with all Gzip-format signature file (eg. *.sig.gz).",
        required=True,
    )
    parser.add_argument(
        "--ksize",
        type=int,
        help="Size of kmers in sketch since Zipfiles can contain multiple k-sizes.",
        required=True,
    )
    parser.add_argument(
        "--num_threads",
        type=_positive_int,
        help="Number of threads to use for parallelization.",
        required=False,
        default=16,
    )
    parser.add_argument(
        "--ani_thresh",
        type=_ani_threshold,
        help="mutation cutoff for species equivalence. Organisms with this ANI "
        'or greater between them are considered "equivalent" (greater than 0, at most 1).',
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
    ref_file = str(Path(args.ref_file).absolute())
    ksize = args.ksize
    num_threads = args.num_threads
    ani_thresh = args.ani_thresh
    prefix = args.prefix
    outdir = str(Path(args.outdir).absolute())
    force = args.force

    # make sure reference database file exist and valid
    logger.info("Checking reference database file")
    if os.path.splitext(ref_file)[1] != ".zip":
        raise ValueError(
            f"Reference database file {ref_file} is not a zip file. Please a Sourmash signature database file with Zipfile format."
        )
    utils.check_file_existence(
        str(Path(ref_file).absolute()),
        f"Reference database zip file {ref_file} does not exist.",
    )

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
                f"Temporary directory {path_to_temp_dir} already exists. Removing it."
            )
            shutil.rmtree(path_to_temp_dir)
    os.makedirs(path_to_temp_dir, exist_ok=True)

    # unzip the sourmash signature file to the temporary directory
    logger.info("Unzipping the sourmash signature file to the temporary directory")
    with zipfile.ZipFile(ref_file, "r") as sourmash_db:
        sourmash_db.extractall(path_to_temp_dir)
    # sourmash stores identical sketches as signatures/<md5>.sig.gz, <md5>.sig.gz_0, <md5>.sig.gz_1, ...
    all_gz_files = [f for f in glob.glob(f"{path_to_temp_dir}/signatures/*") if utils.is_gzipped_sig_file(f)]
    
    # decompress all signature files
    logger.info(f"Decompressing {len(all_gz_files)} .sig.gz files using {num_threads} threads.")
    utils.decompress_all_sig_files(all_gz_files, num_threads)

    # Extract signature information
    logger.info("Extracting signature information")
    sig_info_dict = utils.collect_signature_info(num_threads, ksize, path_to_temp_dir)
    if not sig_info_dict:
        raise ValueError(f"{ref_file} contains no non-empty sketches with k-mer size {ksize}.")
    # check if all signatures have the same ksize and scaled
    logger.info("Checking if all signatures have the same scaled")
    scale_set = set([value[-2] for value in sig_info_dict.values()])
    if len(scale_set) != 1:
        raise ValueError(
            "Not all signatures have the same scaled. Please check your input."
        )
    scale = scale_set.pop()

    # Find the close related genomes with ANI > ani_thresh from the reference database, then remove them, and generate a dataframe with the selected genomes
    logger.info(
        "Finding the closely related genomes with ANI > ani_thresh from the reference database, then remove them, and generate a dataframe with the selected genomes."
    )
    manifest_df = utils.run_yacht_train_core(
        num_threads, ani_thresh, ksize, path_to_temp_dir, sig_info_dict
    )

    # write out the manifest file
    logger.info("Writing out the manifest file")
    manifest_file_path = os.path.join(outdir, f"{prefix}_processed_manifest.tsv")
    manifest_df.to_csv(manifest_file_path, sep="\t", index=None)

    # build the binary index of the selected genomes that `yacht run` uses
    logger.info("Building the YACHT index of the selected genomes")
    index_file_path = os.path.join(outdir, f"{prefix}{INDEX_SUFFIX}")
    utils.build_yacht_index(manifest_df, path_to_temp_dir, ksize, scale, index_file_path)

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
            "index_file_path": index_file_path,
            "yacht_version": YACHT_VERSION,
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
