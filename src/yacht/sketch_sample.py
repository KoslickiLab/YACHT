#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import tempfile
from loguru import logger

# Import global variables

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)



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
    subprocess.run(cmd, shell=True, check=True)
    try:
        logger.info(f"Starting sketching a single-end FASTA/Q file: {infile}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success("Successfully sketched!!")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {infile}: {e}")


def sketch_paired_end(infile1, infile2, kmer, scaled, outfile):
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
        for infile in [infile1, infile2]:
            with open(infile, "r") as f:
                temp_file.write(f.read())
        temp_file_path = temp_file.name

    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund -o {outfile} {temp_file_path}"
    try:
        logger.info(f"Starting sketching paired-end FASTA/Q files: {infile1} {infile2}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success("Successfully sketched!!")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {infile1} {infile2}: {e}")
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
