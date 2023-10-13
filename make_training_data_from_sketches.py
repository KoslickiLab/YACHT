#!/usr/bin/env python
import os, sys
import sourmash
import argparse
import zipfile
from pathlib import Path
import pandas as pd
import srcs.utils as utils
from loguru import logger
import json
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Location of the Sourmash signature database file. '
                                           'This is expected to be in Zipfile format (eg. *.zip)'
                                           'that contains a manifest "SOURMASH-MANIFEST.csv" and a folder "signatures"'
                                           'with all Gzip-format signature file (eg. *.sig.gz) ', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers in sketch since Zipfiles '
                                                  'can contain multiple k-mer sizes', required=True)
    parser.add_argument('--num_threads', type=int, help='Number of threads to use for parallelization.', required=False, default=16)
    parser.add_argument('--ani_thresh', type=float, help='mutation cutoff for species equivalence.',
                        required=False, default=0.95)
    parser.add_argument('--prefix', help='Prefix for this experiment.', required=False, default='yacht')
    parser.add_argument('--outdir', type=str, help='path to output directory', required=False, default=os.getcwd())
    args = parser.parse_args()

    # get the arguments
    ref_file = str(Path(args.ref_file).absolute())
    ksize = args.ksize
    num_threads = args.num_threads
    ani_thresh = args.ani_thresh
    prefix = args.prefix
    outdir = str(Path(args.outdir).absolute())

    # make sure reference database file exist and valid
    logger.info("Checking reference database file")
    if os.path.splitext(ref_file)[1] != '.zip':
        raise ValueError(f"Reference database file {ref_file} is not a zip file. Please a Sourmash signature database file with Zipfile format.")
    utils.check_file_existence(str(Path(ref_file).absolute()), f'Reference database zip file {ref_file} does not exist.')

    # Create a temporary directory with time info as label
    logger.info("Creating a temporary directory")
    path_to_temp_dir = os.path.join(outdir, prefix+'_intermediate_files')
    if os.path.exists(path_to_temp_dir):
        raise ValueError(f"Temporary directory {path_to_temp_dir} already exists. Please remove it or given a new prefix name using parameter '--prefix'.")
    os.makedirs(path_to_temp_dir)
    
    # unzip the sourmash signature file to the temporary directory
    logger.info("Unzipping the sourmash signature file to the temporary directory")
    with zipfile.ZipFile(ref_file, 'r') as sourmash_db:
        sourmash_db.extractall(path_to_temp_dir)

    # Extract signature information
    logger.info("Extracting signature information")
    sig_info_dict = utils.collect_signature_info(num_threads, ksize, path_to_temp_dir)

    # Find the close related genomes with ANI > ani_thresh from the reference database
    logger.info("Find the close related genomes with ANI > ani_thresh from the reference database")
    sig_same_genoms_dict = utils.run_multisearch(num_threads, ani_thresh, ksize, path_to_temp_dir)

    # remove the close related organisms: any organisms with ANI > ani_thresh
    # pick only the one with largest number of unique kmers from all the close related organisms
    logger.info("Removing the close related organisms with ANI > ani_thresh")
    rep_remove_dict, manifest_df = utils.remove_corr_organisms_from_ref(sig_info_dict, sig_same_genoms_dict)

    # write out the manifest file
    logger.info("Writing out the manifest file")
    manifest_file_path = os.path.join(outdir, f'{prefix}_processed_manifest.tsv')
    manifest_df.to_csv(manifest_file_path, sep='\t', index=None)

    # write out a mapping dataframe from representative organism to the close related organisms
    logger.info("Writing out a mapping dataframe from representative organism to the close related organisms")
    rep_remove_df = pd.DataFrame([(rep_org, ','.join(corr_org_list)) for rep_org, corr_org_list in rep_remove_dict.items()])
    rep_remove_df.columns = ['rep_org', 'corr_orgs']
    rep_remove_df_path = os.path.join(outdir, f'{prefix}_rep_to_corr_orgas_mapping.tsv')
    rep_remove_df.to_csv(rep_remove_df_path, sep='\t', index=None)

    # save the config file
    logger.info("Saving the config file")
    json.dump({'manifest_file_path': manifest_file_path,
               'rep_remove_df_path': rep_remove_df_path,
               'pathogen_detection_intermediate_files_dir': path_to_temp_dir,
               'ksize': ksize,
               'ani_thresh': ani_thresh}, open(f'{prefix}_config.json', 'w'), indent=4)
