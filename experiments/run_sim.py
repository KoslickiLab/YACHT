import numpy as np
import os
import sys
SCRIPT_DIR = os.path.abspath(__file__)
sys.path.append(os.path.dirname(SCRIPT_DIR))
import argparse
from Bio import SeqIO
import screed
import pandas as pd
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR))))
import KEGG_sketching_annotation.src.HelperFuncs as KSA


def fasta_to_str(filename):
    """
    Converts a fasta file to a single sequence
    :param filename: location of fasta file
    :return: string consisting of all sequences in file concatenated
    """
    seqs = SeqIO.parse(filename, 'fasta')
    full_seq = ''
    for seq in seqs:
        full_seq += str(seq.seq)
    return full_seq

def format_db_file(genomes_folder, db_filename):
    """
    Rewrites database to bbmap-readable file
    :param genomes_folder: 
    :param db_filename: destination of output file
    :files created: single .fasta file at location db_filename
    :return: None
    """
    genome_dirs = os.listdir(genomes_folder)
    with open(db_filename, "w") as dbf:
        for genome in genome_dirs:
            files = os.listdir(genomes_folder + '/' +genome)
            for filename in files:
                if filename.endswith('_genomic.fna') and not filename.endswith('_from_genomic.fna'):
                    genome_seq = fasta_to_str(genomes_folder + '/' + genome +'/' + filename)
                    dbf.write('>'+filename+'\n')
                    dbf.write(genome_seq+'\n')
                    
def parse_sim(mg_file, count_file):
    """
    Creates csv file with index given by genome fna filenames and one column corresponding to count of reads in sample
    :param mg_file: location of bbmap-produced simulated metagenome file
    :param count_file: location of output .csv file
    :files created: single .csv file at location count_file
    :return: None
    """
    sequence_headers = []
    with screed.open(mg_file) as seqfile:
        for read in seqfile:
            sequence_headers.append(read.name)
            
    simulation_data = {}
    for header in sequence_headers:
        header_split = header.split(".")
        genome_id = header_split[1][1:]+'.'+header_split[2]+'.fna'
        if genome_id not in simulation_data.keys():
            simulation_data[genome_id] = 1
        else:
            simulation_data[genome_id] += 1
    
    count_df = pd.DataFrame(simulation_data.values(), index = simulation_data.keys())
    count_df.to_csv(count_file)
                    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script converts a collection of genome folders into a bbmap-readable database, then uses bbmap to simulate a metagenome.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genomes_folder', help='Location of genome subfolders', required=True)
    parser.add_argument('--out_folder', help='Destination for output files', required=True)
    parser.add_argument('--num_reads', type=int, default=10**6, help='Number of reads to simulate', required=False)
    parser.add_argument('--num_orgs', type=int, default=250, help='Number of organisms to include in the simulation (approximate)', required=False)
    parser.add_argument('--noisy', action='store_true', help='Include this flag to add noise to the simulation', required=False)
    parser.add_argument('--uniform', action='store_true', help='if you want the simulation to be uniform instead '
                                                                   'of exponential', required=False)
    args = parser.parse_args()
    
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)
    
    db_filename = args.out_folder + '/formatted_db.fasta'
    mg_filename = args.out_folder + '/simulated_mg.fq'
    count_filename = args.out_folder + '/simulation_counts.csv'
    if not os.path.exists(db_filename):
        format_db_file(args.genomes_folder, db_filename)
    KSA.run_simulation(db_filename, mg_filename, args.num_reads, num_orgs=args.num_orgs, len_reads=150,
                       noisy=args.noisy, uniform=args.uniform)
    parse_sim(mg_filename, count_filename)
