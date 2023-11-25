#!/usr/bin/env python
import argparse
import subprocess
import os
from pathlib import Path

from yacht import __version__

def sketch_single_file(infile, kmer, scaled, outfile):
    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund --singleton {infile} -o {outfile}"
    subprocess.run(cmd, shell=True, check=True)

def sketch_multiple_files(folder_path, kmer, scaled, outfile):
    dataset_file = os.path.join(folder_path, "dataset.csv")
    with open(dataset_file, "w") as f:
        for path in Path(folder_path).glob('**/*.fasta'):
            f.write(f"{path}\n")
    cmd = f"sourmash sketch fromfile {dataset_file} -p dna,k={kmer},scaled={scaled},abund -o {outfile}"
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Sketch genomes using Sourmash.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument("--infile", help="Input file or folder path.", required=True)
    parser.add_argument("--kmer", type=int, help="K-mer size.", default=31)
    parser.add_argument("--scaled", type=int, help="Scaled factor.", default=1000)
    parser.add_argument("--outfile", help="Output file name.", required=True)
    args = parser.parse_args()

    # Check if input is a file or a directory
    if os.path.isfile(args.infile):
        sketch_single_file(args.infile, args.kmer, args.scaled, args.outfile)
    elif os.path.isdir(args.infile):
        sketch_multiple_files(args.infile, args.kmer, args.scaled, args.outfile)
    else:
        raise FileNotFoundError(f"Input path {args.infile} does not exist.")

if __name__ == "__main__":
    main()
