#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
from pathlib import Path
from loguru import logger
# Import global variables

# Configure Loguru logger
logger.remove()
logger.add(
    sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO"
)


def add_arguments(parser):
    parser.add_argument("--infile", help="Input file or folder path.", required=True)
    parser.add_argument("--kmer", type=int, help="K-mer size.", default=31)
    parser.add_argument("--scaled", type=int, help="Scaled factor.", default=1000)
    parser.add_argument("--outfile", help="Output file name.", required=True)


def sketch_single_file(infile, kmer, scaled, outfile):
    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund --singleton -o {outfile} {infile}"
    try:
        logger.info(f"Starting sketching a single file: {infile}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success("Successfully sketched!!")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {infile}: {e}")


def sketch_multiple_files(folder_path, kmer, scaled, outfile):
    dataset_file = os.path.join(folder_path, "dataset.csv")
    file_extensions = [
        "*.fasta",
        "*.fna",
        "*.fas",
        "*.fa",
        "*.fasta.gz",
        "*.fna.gz",
        "*.fas.gz",
        "*.fa.gz",
    ]

    try:
        logger.info(
            f"Preparing dataset file for multiple sequence files in {folder_path}"
        )
        with open(dataset_file, "w") as f:
            f.write("name,genome_filename,protein_filename\n")  # Add header for CSV
            for extension in file_extensions:
                for path in Path(folder_path).glob(f"**/{extension}"):
                    # Use the file name (without extension) as the sketch name
                    name = path.name.replace(extension.replace("*", ""), "")
                    f.write(
                        f"{name},{str(path.absolute())},\n"
                    )  # Write the name and the full file path

        cmd = f"sourmash sketch fromfile {dataset_file} -p dna,k={kmer},scaled={scaled},abund -o {outfile} --force-output-already-exists"
        logger.info(f"Starting sketching multiple sequence files in: {folder_path}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched files in: {folder_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching files in {folder_path}: {e}")


def main(args):
    try:
        if os.path.isfile(args.infile):
            sketch_single_file(args.infile, args.kmer, args.scaled, args.outfile)
        elif os.path.isdir(args.infile):
            sketch_multiple_files(args.infile, args.kmer, args.scaled, args.outfile)
        else:
            raise FileNotFoundError(f"Input path {args.infile} does not exist.")
    except FileNotFoundError as e:
        logger.error(e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sketch genomes using Sourmash.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)
