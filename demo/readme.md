### A toy sample for YACHT demonstration
How to get the 2 sketch files here?
Randomly download 20 genomes from GTDB into the current directory.

```
# build ref for 20 genomes
sourmash sketch dna -f -p k=31,scaled=1000,abund -o ref.sig.zip *.fna.gz

# simulate metagenomics for the 1st 10 genomes and sketch
cat $(ls *fna.gz | head) > merged_10.fna.gz
randomreads.sh ref=merged_10.fna.gz out=BBMap_simulated_meta.fq coverage=3 len=100 metagenome
sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip BBMap_simulated_meta.fq
```





