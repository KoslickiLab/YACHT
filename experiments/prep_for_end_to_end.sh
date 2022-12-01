python ../../KEGG_sketching_annotation/scripts/create_genome_ref_db.py ./ref_genomes_3/reference_genomes 1045_database 1045
# run the sim once
python run_sim.py --genomes_folder ref_genomes_3/reference_genomes --out_folder ${simsFolder} --num_reads 1000 --num_orgs 1000
# then copy this to the `experiments` dir
cp sims-*/formatted_db.fasta ..
# sketch the formatted_db.fasta
sourmash sketch dna -f -p k=31,scaled=1000,abund -o formatted_db.sig --singleton formatted_db.fasta
# and create a manifest for it
sourmash sig manifest formatted_db.sig -o MANIFEST.csv
