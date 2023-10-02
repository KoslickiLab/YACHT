#!/usr/bin/env python
import sys
import sourmash
import argparse
from pathlib import Path
from scipy.sparse import save_npz
import srcs.utils as utils
from loguru import logger
import json
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script converts a collection of signature files into a reference database matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_file', help='Location of the Sourmash signature file. '
                                           'This is expected to be in Zipfile format (eg. *.sig.zip)', required=True)
    parser.add_argument('--ksize', type=int, help='Size of kmers in sketch since Zipfiles '
                                                  'can contain multiple k-mer sizes', required=True)
    parser.add_argument('--ani_thresh', type=float, help='mutation cutoff for species equivalence.',
                        required=False, default=0.95)
    parser.add_argument('--out_prefix', help='Location and prefix for output files.', required=True)
    args = parser.parse_args()

    # get the arguments
    ref_file = args.ref_file
    ksize = args.ksize
    ani_thresh = args.ani_thresh
    out_prefix = args.out_prefix

    # load the signatures
    logger.info(f"Loading signatures from {ref_file}")
    signatures = sourmash.load_file_as_signatures(ref_file)
    signature_count = utils.count_files_in_zip(ref_file) - 1

    # DONE: do signature size checking, coverting to sourmash list and generate reference matrix at the same time
    # check that all signatures have the same ksize as the one provided
    # signatures_mismatch_ksize return False (if all signatures have the same kmer size)
    # or True (the first signature with a different kmer size)
    # convert signatures to reference matrix (rows are hashes/kmers, columns are organisms)
    logger.info("Converting signatures to reference matrix")
    signatures, ref_matrix, hashes, is_mismatch = utils.signatures_to_ref_matrix(signatures, ksize, signature_count)
    if is_mismatch:
        raise ValueError(f"Not all signatures from sourmash signature file {ref_file} have the given ksize {ksize}")

    # remove 'same' organisms: any organisms with ANI > ani_thresh are considered the same organism
    logger.info("Removing 'same' organisms with ANI > ani_thresh")
    processed_ref_matrix, uncorr_org_idx = utils.get_uncorr_ref(ref_matrix, ksize, ani_thresh)
    save_npz(f'{out_prefix}_ref_matrix_processed.npz', processed_ref_matrix)

    # write out hash-to-row-indices file
    logger.info("Writing out hash-to-row-indices file")
    utils.write_hashes(f'{out_prefix}_hash_to_col_idx.pkl', hashes)

    # write out organism manifest (original index, processed index, num unique kmers, num total kmers, scale factor)
    logger.info("Writing out organism manifest")
    utils.write_processed_indices(f'{out_prefix}_processed_org_idx.csv', signatures, uncorr_org_idx)

    # save the k-mer size and ani threshold to a json file
    logger.info("Saving k-mer size and ani threshold to json file")
    json.dump({'reference_matrix_path': str(Path(f'{out_prefix}_ref_matrix_processed.npz').resolve()),
               'hash_to_idx_path': str(Path(f'{out_prefix}_hash_to_col_idx.pkl').resolve()),
               'processed_org_file_path': str(Path(f'{out_prefix}_processed_org_idx.csv').resolve()),
               'ksize': ksize,
               'ani_thresh': ani_thresh}, open(f'{out_prefix}_config.json', 'w'), indent=4)
