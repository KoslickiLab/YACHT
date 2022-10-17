import csv

def load_hashes(filename):
    with open(filename, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        hashes = {int(rows[0]):int(rows[1]) for rows in reader}
    return hashes

def load_processed_organisms(filename):
    with open(filename, mode='r') as infile:
        next(infile)
        reader = csv.reader(infile)
        orgs = [rows[0] for rows in reader]
    return orgs