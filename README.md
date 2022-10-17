# Estimating_Unknowns
Creating a reference dictionary matrix (ref_matrix.py):

python ref_matrix.py --ref_file '../ForSteve/ref_gtdb-rs207.genomic-reps.dna.k31.zip' --out_prefix 'test2_' --N 20


Computing relative abundance of organisms (recover_abundance.py):

python recover_abundance.py --ref_file 'test2_ref_matrix_processed.npz' --sample_file '../ForSteve/sample.sig' --hash_file 'test2_hash_to_col_idx.csv' --org_file 'test2_processed_org_idx.csv' --w 0.01 --outfile 'test2_recovered_abundance.csv'
