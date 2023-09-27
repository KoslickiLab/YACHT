These signatures are obtained from /data/jzr5814/repositories/YACHT/tests/testdata directory. I don't know what is in these signatures but we can take a look at what's inthese signatures with the following commands.


```bash
sourmash signature fileinfo sample.sig 
```

== This is sourmash version 4.8.3. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

** loading from 'sample.sig'
path filetype: MultiIndex
location: sample.sig
is database? no
has manifest? yes
num signatures: 1
** examining manifest...
total hashes: 49821
summary of sketches:
   1 sketches with DNA, k=31, scaled=1000, abund      49821 total hashes

```bash
sourmash signature fileinfo 20_genomes_sketches.zip 
```

== This is sourmash version 4.8.3. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

** loading from '20_genomes_sketches.zip'
path filetype: ZipFileLinearIndex
location: /data/jzr5814/repositories/YACHT/use_case_examples/example_1/20_genomes_sketches.zip
is database? yes
has manifest? yes
num signatures: 20
** examining manifest...
total hashes: 63898
summary of sketches:
   20 sketches with DNA, k=31, scaled=1000, abund     63898 total hashes

```bash
python ../../make_training_data_from_sketches.py --ref_file '20_genomes_sketches.zip' --ksize 31 --out_prefix '20_genomes' --ani_thresh 0.95
```

023-09-27 08:54:59 - INFO - Loading signatures from 20_genomes_sketches.zip
2023-09-27 08:54:59 - INFO - Converting signatures to reference matrix
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 855.31it/s]
2023-09-27 08:54:59 - INFO - Removing 'same' organisms with ANI > ani_thresh
2023-09-27 08:54:59 - INFO - Writing out hash-to-row-indices file
2023-09-27 08:54:59 - INFO - Writing out organism manifest
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:00<00:00, 570.85it/s]
2023-09-27 08:54:59 - INFO - Saving k-mer size and ani threshold to json file

