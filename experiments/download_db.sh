python KEGG_sketching_annotation/scripts/get_reference_genomes.py -n 100 -u -s ref_genomes_3 -S 2 -c sew347@psu.edu

python KEGG_sketching_annotation/scripts/create_genome_ref_db.py --in_dr ref_genomes_3/reference_genomes --db_name ref_genomes/test_database --num_genomes 100