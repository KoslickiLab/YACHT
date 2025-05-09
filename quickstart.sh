NUM_THREADS=1 # Adjust based on your machine's capabilities

cd demo # the 'demo' folder can be downloaded via command 'yacht download demo' if it doesn't exist

# # build k-mer sketches for the query sample and ref genomes
# yacht sketch sample --infile ./query_data/query_data.fq --kmer 31 --scaled 1000 --outfile sample.sig.zip
# yacht sketch ref --infile ./ref_genomes --kmer 31 --scaled 1000 --outfile ref.sig.zip

# # preprocess the reference genomes (training step)
yacht train --ref_file ref.sig.zip --ksize 31 --num_threads ${NUM_THREADS} --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./ --force

# run YACHT algorithm to check the presence of reference genomes in the query sample (inference step)
yacht run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads ${NUM_THREADS} --min_coverage_list 1 0.6 0.2 0.1 --out ./result.xlsx

# # convert result to CAMI profile format (Optional)
# yacht convert --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./


