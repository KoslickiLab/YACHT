# YACHT

YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity. 

**_PLEASE NOTE: the proof-of-concept code is contained [here](https://github.com/KoslickiLab/YACHT-Proof-of-concept). This repo is a work-in-progress for the production level YACHT implementation._**

## Installation
```bash
conda install -c conda-forge -c bioconda -c anaconda sourmash=4.5.0 cvxpy scipy numpy pandas scikit-learn
```

## Creating a reference dictionary matrix (`ref_matrix.py`):
```bash 
python ref_matrix.py --ref_file '../ForSteve/ref_gtdb-rs207.genomic-reps.dna.k31.zip' --out_prefix 'test2_' --N 20
```

## Computing relative abundance of organisms (`recover_abundance.py`):
```bash
python recover_abundance.py --ref_file 'test2_ref_matrix_processed.npz' --sample_file '../ForSteve/sample.sig' --hash_file 'test2_hash_to_col_idx.csv' --org_file 'test2_processed_org_idx.csv' --w 0.01 --outfile 'test2_recovered_abundance.csv'
```

## Basic workflow
1. run ```python ref_matrix.py --ref_file 'tests/testdata/20_genomes_sketches.zip' --out_prefix 'tests/unittest_'``` . This 
should generate 4 files in the tests folder.
2. run ```python recover_abundance.py --ref_file 'tests/unittest_ref_matrix_processed.npz' --sample_file 
   'tests/testdata/sample.sig' --hash_file 'tests/unittest_hash_to_col_idx.csv' --org_file 'tests/unittest_processed_org_idx.csv' --w 0.01 --outfile 'tests/unittest_recovered_abundance.csv'``` . Should create a file `tests/unittest_recovered_abundance.csv` which should be all zeros.
3. run the same command as above, but with `--w 0.0001`. Should overwrite `tests/unittest_recovered_abundance.csv` with a 
   6 in the 19th row
