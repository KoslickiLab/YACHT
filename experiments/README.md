# What does what

## `prep_for_end_to_end.sh`
Get things ready for running `end_to_end.sh`.

## `end_to_end.sh` and `end_to_end_51.sh`
Run the end to end pipeline for k=31 or k=51. In short, it will make a  folder `sims-<num>` where `<num>` is a time 
stamp. Results are in that folder in a file `results.txt` and appends the results to `../all_results.txt`. Don't run too 
many of them in parallel, since they make some pretty massive intermediate files. Important parameters include

1. `numReads`: number of reads to simulate
2. `--num_orgs`: number of organisms to include in the simulation
2. `N`: number of organisms you would like to have marked as _potentially unknown_. The true number of unknowns is 
calculated with `calculate_unknown_percent.py`


## `single_genome.sh`
This will run a "simulation" where the reference is a single real genome, and the query is a single real genome with 
a known distance to the reference. This relies on the `one_genome_sample.py` file that gets everthing needed by the 
guts in `end_to_end.sh`. Results are stored in the file `results.txt` and _not_ appended to `all_results.txt`. 
Important parameters include:

1. `--out_dir` where you want this to take place (use something like `sims-$(date +"%s")` if you want a date stamp). 
   Defaults to `one_genome_sample`
2. `--bucket_number` ANI values are split into buckets of width 0.01 ANI, with the first bucket being the pairs of 
   genomes that have [1, 0.99) ANI. This parameter specifies which bucket to use. Currently set to 0 (the highest 
   ANI bucket). Looks like the buckets become empty after the 26th one: [0.75, 0.74).
3. `--pair_number` Within each bucket, there are some variable number of pairs of genomes. This parameter specifies 
   which pair to use. Currently set to 0 (the first pair in the bucket). Looks like there is a minimum of 24 pairs 
   in each bucket (though some have a ton more). Change this if you want to see what happens if you select a 
   different pair of genomes with the same ANI.

