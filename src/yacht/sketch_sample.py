#!/usr/bin/env python
import argparse
import subprocess
import os
import re
import sys
import tempfile
from importlib import metadata
from loguru import logger

# Import global variables

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)

# `sourmash scripts singlesketch` (sourmash_plugin_branchwater >= 0.9.13, which accepts several input files) sketches
# with sourmash's own Rust code: the same hashes and abundances as `sourmash sketch dna`, without its per-read Python
# overhead and without concatenating paired-end files first. Older plugin versions and other output formats use
# `sourmash sketch dna`.
SINGLESKETCH_MIN_VERSION = (0, 9, 13)
SINGLESKETCH_OUTPUT_SUFFIXES = (".zip", ".sig", ".sig.gz")


def use_singlesketch(outfile):
    if not str(outfile).endswith(SINGLESKETCH_OUTPUT_SUFFIXES):
        return False
    try:
        version = re.match(r"(\d+)\.(\d+)\.(\d+)", metadata.version("sourmash_plugin_branchwater"))
    except metadata.PackageNotFoundError:
        return False
    return version is not None and tuple(int(x) for x in version.groups()) >= SINGLESKETCH_MIN_VERSION


def singlesketch(infiles, kmer, scaled, outfile):
    """
    Sketch infiles into one signature with `sourmash scripts singlesketch`.
    :return: False if singlesketch rejected its command line (the caller then uses `sourmash sketch dna`)
    """
    # run the sourmash of this Python environment, whose plugin version use_singlesketch checked
    cmd = [sys.executable, "-m", "sourmash", "scripts", "singlesketch", "-p", f"k={kmer},scaled={scaled},abund"]
    cmd += ["-o", str(outfile)]
    cmd += [str(infile) for infile in infiles]
    # singlesketch runs in a single thread: keep the libraries it loads from starting idle thread pools
    env = dict(os.environ)
    for var in ("RAYON_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env.setdefault(var, "1")
    result = subprocess.run(cmd, env=env)
    if result.returncode == 2:  # argparse usage error, e.g. a plugin release with a different command line
        logger.warning("`sourmash scripts singlesketch` rejected its arguments, using `sourmash sketch dna` instead.")
        return False
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return True


def add_arguments(parser):
    parser.add_argument(
        "--infile",
        nargs="+",
        help="Input FASTA/Q file(s). For paired-end reads, provide two files.",
        required=True,
    )
    parser.add_argument("--kmer", type=int, help="K-mer size.", default=31)
    parser.add_argument("--scaled", type=int, help="Scaled factor.", default=1000)
    parser.add_argument("--outfile", help="Output file name.", required=True)


def sketch_single_end(infile, kmer, scaled, outfile):
    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund -o {outfile} {infile}"
    try:
        logger.info(f"Starting sketching a single-end FASTA/Q file: {infile}")
        if not (use_singlesketch(outfile) and singlesketch([infile], kmer, scaled, outfile)):
            subprocess.run(cmd, shell=True, check=True)
        logger.success("Successfully sketched!!")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {infile}: {e}")
        raise


def sketch_paired_end(infile1, infile2, kmer, scaled, outfile):
    try:
        logger.info(f"Starting sketching paired-end FASTA/Q files: {infile1} {infile2}")
        if not (use_singlesketch(outfile) and singlesketch([infile1, infile2], kmer, scaled, outfile)):
            sketch_concatenated_pair(infile1, infile2, kmer, scaled, outfile)
        logger.success("Successfully sketched!!")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {infile1} {infile2}: {e}")
        raise


def sketch_concatenated_pair(infile1, infile2, kmer, scaled, outfile):
    # `sourmash sketch dna` makes one signature per input file, so concatenate the pair first
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        for infile in [infile1, infile2]:
            with open(infile, "r") as f:
                temp_file.write(f.read())
        temp_file_path = temp_file.name

    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund -o {outfile} {temp_file_path}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    finally:
        os.remove(temp_file_path)


def main(args):
    if len(args.infile) == 1:
        sketch_single_end(args.infile[0], args.kmer, args.scaled, args.outfile)
    elif len(args.infile) == 2:
        sketch_paired_end(
            args.infile[0], args.infile[1], args.kmer, args.scaled, args.outfile
        )
    else:
        raise ValueError(
            "Please provide either one file for single-end reads or two files for paired-end reads."
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sketch single-end or paired-end reads using Sourmash.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
