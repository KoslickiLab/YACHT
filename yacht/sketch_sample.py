#!/usr/bin/env python
import argparse
import subprocess
import os
from pathlib import Path
import tempfile

from yacht import __version__

def sketch_single_end(infile, kmer, scaled, outfile):
    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund -o {outfile} {infile}"
    subprocess.run(cmd, shell=True, check=True)

def sketch_paired_end(infile1, infile2, kmer, scaled, outfile):
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        for infile in [infile1, infile2]:
            with open(infile, 'r') as f:
                temp_file.write(f.read())
        temp_file_path = temp_file.name

    cmd = f"sourmash sketch dna -f -p k={kmer},scaled={scaled},abund -o {outfile} {temp_file_path}"
    subprocess.run(cmd, shell=True, check=True)
    os.remove(temp_file_path)

def main():
    parser = argparse.ArgumentParser(description="Sketch single-end or paired-end reads using Sourmash.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument("--infile", nargs='+', help="Input FASTA/Q file(s). For paired-end reads, provide two files.", required=True)
    parser.add_argument("--kmer", type=int, help="K-mer size.", default=31)
    parser.add_argument("--scaled", type=int, help="Scaled factor.", default=1000)
    parser.add_argument("--outfile", help="Output file name.", required=True)
    args = parser.parse_args()

    if len(args.infile) == 1:
        sketch_single_end(args.infile[0], args.kmer, args.scaled, args.outfile)
    elif len(args.infile) == 2:
        sketch_paired_end(args.infile[0], args.infile[1], args.kmer, args.scaled, args.outfile)
    else:
        raise ValueError("Please provide either one file for single-end reads or two files for paired-end reads.")

if __name__ == "__main__":
    main()

