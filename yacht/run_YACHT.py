#!/usr/bin/env python
import os
import sys
import pandas as pd
from pathlib import Path
from . import hypothesis_recovery_src as hr
import argparse
from . import utils
import json
import warnings
import zipfile
from loguru import logger

warnings.filterwarnings("ignore")

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)


def add_arguments(parser):
    parser.add_argument(
        "--version", action="version", version=f"YACHT {utils.__version__}"
    )
    parser.add_argument(
        "--json",
        type=str,
        help="Path to a json file generated by make_training_data_from_sketches.py.",
        required=True,
    )
    parser.add_argument(
        "--sample_file", help="Metagenomic sample in .sig.zip format", required=True
    )
    parser.add_argument(
        "--significance",
        type=float,
        help="Minimum probability of individual true negative.",
        required=False,
        default=0.99,
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        help="Number of threads to use for parallelization.",
        required=False,
        default=16,
    )
    parser.add_argument(
        "--keep_raw", action="store_true", help="Keep raw results in output file."
    )
    parser.add_argument(
        "--show_all",
        action="store_true",
        help="Show all organisms (no matter if present) in output file.",
    )
    parser.add_argument(
        "--min_coverage_list",
        nargs="+",
        type=float,
        help="A list of percentages of unique k-mers covered by reads in the sample. "
        "Each value should be between 0 and 1, with 0 being the most sensitive (and least "
        "precise) and 1 being the most precise (and least sensitive).",
        required=False,
        default=[1, 0.5, 0.1, 0.05, 0.01],
    )
    parser.add_argument(
        "--out",
        type=str,
        help="path to output excel file",
        required=False,
        default=os.path.join(os.getcwd(), "result.xlsx"),
    )


def main(args):
    json_file_path = str(Path(args.json).absolute())  # path to json file
    sample_file = str(Path(args.sample_file).absolute())  # location of sample.sig file
    significance = args.significance  # Minimum probability of individual true negative.
    num_threads = args.num_threads  # Number of threads to use for parallelization.
    keep_raw = args.keep_raw  # Keep raw results in output file.
    show_all = args.show_all # Show all organisms (no matter if present) in output file.
    min_coverage_list = args.min_coverage_list # a list of percentages of unique k-mers covered by reads in the sample.
    out = str(Path(args.out).absolute())  # full path to output excel file
    outdir = os.path.dirname(out)  # path to output directory
    out_filename = os.path.basename(out)  # output filename

    # check if the output filename is valid
    if os.path.splitext(out_filename)[1] != ".xlsx":
        raise ValueError(
            f"Output filename {out} is not a valid excel file. Please use .xlsx as the extension."
        )

    # check if the json file exists
    utils.check_file_existence(
        json_file_path,
        f"Config file {json_file_path} does not exist. "
        f"Please run make_training_data_from_sketches.py first.",
    )
    # load the config file, ksize, and ani_thresh
    config = json.load(open(json_file_path, "r"))
    manifest_file_path = config["manifest_file_path"]
    path_to_genome_temp_dir = config["intermediate_files_dir"]
    scale = config["scale"]
    ksize = config["ksize"]
    ani_thresh = config["ani_thresh"]

    # Make sure the output can be written to
    if not os.access(outdir, os.W_OK):
        # give error message, and exit with error status
        print(f"Cannot write to the location: {outdir}.\n")
        print("Please check if this location exists, and that you have the permission to write to this location. Exiting..\n")
        sys.exit(1)

    # check if min_coverage is between 0 and 1
    for x in min_coverage_list:
        if not (0 <= x <= 1):
            raise ValueError(
                f"One of values in the min_coverage_list you provided {x} is not between 0 and 1. Please check your input."
            )

    # make sure all these files exist
    utils.check_file_existence(
        manifest_file_path,
        f"The manifest file {manifest_file_path} "
        f"does not exist. Please check if you are using the correct json file as input.",
    )

    # load the training data
    logger.info("Loading the manifest file generated from the training data.")
    manifest = pd.read_csv(manifest_file_path, sep="\t", header=0)

    # check if there is a manifest in the sample sig file
    with zipfile.ZipFile(sample_file, "r") as zip_file:
        if "SOURMASH-MANIFEST.csv" not in zip_file.namelist():
            raise FileNotFoundError(
                f"The input file {sample_file} appears to be missing a manifest associated with it. "
                f"Try running: sourmash sig merge {sample_file} -o <new signature with the manifest present>. "
                f"And then run YACHT using the output of that command."
            )

    # load sample signature and its signature info
    logger.info("Loading sample signature and its signature info.")
    try:
        sample_sig = utils.load_signature_with_ksize(sample_file, ksize)
    except ValueError:
        raise ValueError(
            f"Expected exactly one signature with ksize {ksize} in {sample_file}, found {len(sample_file)}. "
            f"Likely you will need to do something like: sourmash sig merge {sample_file} -o <new signature with just one sketch in it>."
        )
    sample_sig_info = utils.get_info_from_single_sig(sample_file, ksize)

    # add sample signature info to the manifest
    manifest["num_exclusive_kmers_in_sample_sketch"] = sample_sig_info[3]
    manifest["num_total_kmers_in_sample_sketch"] = utils.get_num_kmers(
        sample_sig_info[2], sample_sig_info[3], sample_sig_info[4], scale=False
    )
    manifest["sample_scale_factor"] = sample_sig_info[4]
    manifest["min_coverage"] = 1.0

    # check that the sample scale factor is the same as the genome scale factor for all organisms
    if scale != sample_sig_info[4]:
        raise ValueError(
            "Sample scale factor does not equal genome scale factor. Please check your input."
        )

    # compute hypothesis recovery
    logger.info("Computing hypothesis recovery.")
    sample_info_set = (sample_file, sample_sig)
    min_coverage_list = list(set(min_coverage_list))
    min_coverage_list.sort(reverse=True)
    has_raw = False
    if 1.0 not in min_coverage_list:
        min_coverage_list = [1.0] + min_coverage_list
    else:
        has_raw = True

    manifest_list = hr.hypothesis_recovery(
        manifest,
        sample_info_set,
        path_to_genome_temp_dir,
        min_coverage_list,
        scale,
        ksize,
        significance,
        ani_thresh,
        num_threads,
    )

    # remove unnecessary columns
    remove_cols = ["md5sum", "sample_scale_factor"]
    temp_manifest_list = []
    for temp_manifest in manifest_list:
        temp_manifest = temp_manifest[
            [col for col in temp_manifest.columns if col not in remove_cols]
        ]
        temp_manifest.rename(
            columns={"genome_scale_factor": "scale_factor"}, inplace=True
        )
        temp_manifest_list += [temp_manifest]
    manifest_list = temp_manifest_list

    # process the results and save them to an Excel file
    logger.info(f"Saving results to {outdir}.")
    # save the results with different min_coverage
    with pd.ExcelWriter(out, engine="openpyxl", mode="w") as writer:
        # save the raw results (i.e., min_coverage=1.0)
        if keep_raw:
            temp_mainifest = manifest_list[0].copy()
            temp_mainifest.rename(
                columns={
                    "acceptance_threshold_with_coverage": "acceptance_threshold_wo_coverage",
                    "actual_confidence_with_coverage": "actual_confidence_wo_coverage",
                    "alt_confidence_mut_rate_with_coverage": "alt_confidence_mut_rate_wo_coverage",
                },
                inplace=True,
            )
            temp_mainifest.to_excel(writer, sheet_name="raw_result", index=False)
        # save the results with different min_coverage given by the user
        if not has_raw:
            min_coverage_list = min_coverage_list[1:]
            manifest_list = manifest_list[1:]

        for min_coverage, temp_mainifest in zip(min_coverage_list, manifest_list):
            if not show_all:
                temp_mainifest = temp_mainifest[temp_mainifest["in_sample_est"] == True]
            temp_mainifest.to_excel(
                writer, sheet_name=f"min_coverage{min_coverage}", index=False
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script estimates the abundance of microorganisms from a "
        "reference database matrix and metagenomic sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
