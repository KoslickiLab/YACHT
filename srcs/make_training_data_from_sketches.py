#!/usr/bin/env python
import os, sys
import argparse
import zipfile
from pathlib import Path
from loguru import logger
import json
import shutil

from srcs import utils

logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

def add_arguments(parser):
    parser.add_argument('--ref_file', help='Location of the Sourmash signature database file. '
                                           'This is expected to be in Zipfile format (eg. *.zip)'
                                           'that contains a manifest "SOURMASH-MANIFEST.csv" and a folder "signatures"'
                                           'with all Gzip-format signature file (eg. *.sig.gz) ', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers in sketch since Zipfiles', required=True)
    parser.add_argument('--num_threads', type=int, help='Number of threads to use for parallelization.', required=False,
                        default=16)
    parser.add_argument('--ani_thresh', type=float, help='mutation cutoff for species equivalence.',
                        required=False, default=0.95)
    parser.add_argument('--prefix', help='Prefix for this experiment.', required=False, default='yacht')
    parser.add_argument('--outdir', type=str, help='path to output directory', required=False, default=os.getcwd())
    parser.add_argument('--force', action='store_true', help='Overwrite the output directory if it exists')


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
    if os.path.splitext(ref_file)[1] != '.zip':
        raise ValueError(
            f"Reference database file {ref_file} is not a zip file. Please a Sourmash signature database file with Zipfile format.")
    utils.check_file_existence(str(Path(ref_file).absolute()),
                               f'Reference database zip file {ref_file} does not exist.')

    # Create a temporary directory with time info as label
    logger.info("Creating a temporary directory")
    path_to_temp_dir = os.path.join(outdir, prefix + '_intermediate_files')
    if os.path.exists(path_to_temp_dir) and not force:
        raise ValueError(
            f"Temporary directory {path_to_temp_dir} already exists. Please remove it or given a new prefix name using parameter '--prefix'.")
    else:
        # remove the temporary directory if it exists
        if os.path.exists(path_to_temp_dir):
            logger.warning(f"Temporary directory {path_to_temp_dir} already exists. Removing it.")
            shutil.rmtree(path_to_temp_dir)
    os.makedirs(path_to_temp_dir, exist_ok=True)

    # unzip the sourmash signature file to the temporary directory
    logger.info("Unzipping the sourmash signature file to the temporary directory")
    with zipfile.ZipFile(ref_file, 'r') as sourmash_db:
        sourmash_db.extractall(path_to_temp_dir)

    # Extract signature information
    logger.info("Extracting signature information")
    sig_info_dict = utils.collect_signature_info(num_threads, ksize, path_to_temp_dir)
    # check if all signatures have the same ksize and scaled
    logger.info("Checking if all signatures have the same scaled")
    scale_set = set([value[-1] for value in sig_info_dict.values()])
    if len(scale_set) != 1:
        raise ValueError(f"Not all signatures have the same scaled. Please check your input.")
    scale = scale_set.pop()

    # Find the close related genomes with ANI > ani_thresh from the reference database
    logger.info("Find the close related genomes with ANI > ani_thresh from the reference database")
    multisearch_result = utils.run_multisearch(num_threads, ani_thresh, ksize, scale, path_to_temp_dir)

    # remove the close related organisms: any organisms with ANI > ani_thresh
    # pick only the one with largest number of unique kmers from all the close related organisms
    logger.info("Removing the close related organisms with ANI > ani_thresh")
    remove_corr_df, manifest_df = utils.remove_corr_organisms_from_ref(sig_info_dict, multisearch_result)

    # write out the manifest file
    logger.info("Writing out the manifest file")
    manifest_file_path = os.path.join(outdir, f'{prefix}_processed_manifest.tsv')
    manifest_df.to_csv(manifest_file_path, sep='\t', index=None)

    # write out a mapping dataframe from representative organism to the close related organisms
    logger.info("Writing out a mapping dataframe from representative organism to the close related organisms")
    if len(remove_corr_df) == 0:
        logger.warning("No close related organisms found.")
        remove_corr_df_indicator = ""
    else:
        remove_corr_df_path = os.path.join(outdir, f'{prefix}_removed_orgs_to_corr_orgas_mapping.tsv')
        remove_corr_df.to_csv(remove_corr_df_path, sep='\t', index=None)
        remove_corr_df_indicator = remove_corr_df_path

    # save the config file
    logger.info("Saving the config file")
    json_file_path = os.path.join(outdir, f'{prefix}_config.json')
    json.dump({'manifest_file_path': manifest_file_path,
               'remove_cor_df_path': remove_corr_df_indicator,
               'intermediate_files_dir': path_to_temp_dir,
               'scale': scale,
               'ksize': ksize,
               'ani_thresh': ani_thresh}, open(json_file_path, 'w'), indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
