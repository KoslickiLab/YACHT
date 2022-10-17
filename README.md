# Estimating_Unknowns

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