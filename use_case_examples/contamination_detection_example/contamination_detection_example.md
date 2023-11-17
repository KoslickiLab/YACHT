# Contamination Detection Example
A research is being conducted on how microbial communities are being shaped among diffrent types of respiratory diseases. Samples were collected from two patience in wihch one patient is M. pneumoniae positive and the other in H. influenzae. To save time and many, only one 96-well tray will be used for both samples. Before downstream analysis can be performed, we want to know if cross contamination between samples occured during the loading of the 96-well tray and we randomly choose wells 11, 23, 64, and 80.

Make sure all bacterial reads needed to create your reference dataset also known as a training dataset are available.
```bash
bash 1_before_starting.sh
```
```bash
bash 2_before_starting.sh
```

### Sketch your training dataset and sample to your preference.

#### Using k=31
Note: training and sample datasets are required to have the same ksize. Please note that since we are sketching from a list of genomes. We can use the following sourmash sketch command:
```bash
sourmash sketch fromfile genome_list.csv -p dna,k=31,scaled=1000,abund -o training_database.k31.sig.zip
```

Sketch the negative control reads from well 11
```bash
sourmash sketch dna negative_control_well_11.fasta -p k=31,scaled=1000,abund -o negative_control_well_11.k31.sig.zip
```

Sketch the positive control from well 23
```bash
sourmash sketch dna positive_control_well_23.fasta -p k=31,scaled=1000,abund -o positive_control_well_23.k31.sig.zip
```

Sketch the positive control from well 64
```bash
sourmash sketch dna positive_control_well_64.fasta -p k=31,scaled=1000,abund -o positive_control_well_64.k31.sig.zip
```

Sketch the sample from well 80
```bash
sourmash sketch dna sample_well_80.fasta -p k=31,scaled=1000,abund -o sample_well_80.k31.sig.zip
```

### Make training data for k=31
```bash
python ../../make_training_data_from_sketches.py --ref_file training_database.k31.sig.zip --ksize 31 --ani_thresh 0.95 --out_prefix 'training_database.k31'
```

### Identify whether the patient has a infection and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json 'training_database.k31_config.json' --sample_file 'negative_control_well_11.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'negative_control_well_11_k31_result.xlsx' --outdir './'
```

```bash
python ../../run_YACHT.py --json 'training_database.k31_config.json' --sample_file 'positive_control_well_23.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'positive_control_well_23_k31_result.xlsx' --outdir './'
```

```bash
python ../../run_YACHT.py --json 'training_database.k31_config.json' --sample_file 'positive_control_well_64.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'positive_control_well_64_k31_result.xlsx' --outdir './'
```

```bash
python ../../run_YACHT.py --json 'training_database.k31_config.json' --sample_file 'sample_well_80.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'sample_well_80_k31_result.xlsx' --outdir './'
```

### Results
Using a ksize of 31 at ANI 0.95, YACHT finds XYZ

## Let's decrease ANI to 0.50

### Make training data for k=31
```bash
python ../../make_training_data_from_sketches.py --ref_file training_database.k31.sig.zip --ksize 31 --ani_thresh 0.50 --out_prefix 'training_database.k31_ani0.50'
```

### Pathogen Detection using YACHT
Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json 'training_database.k31_ani0.50_config.json' --sample_file 'negative_control_well_11.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'k31_ani0.50_result.xlsx' --outdir './'
```

Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json 'training_database.k31_ani0.50_config.json' --sample_file 'positive_control_well_23.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'k31_ani0.50_result.xlsx' --outdir './'
```

Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json 'training_database.k31_ani0.50_config.json' --sample_file 'positive_control_well_64.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'k31_ani0.50_result.xlsx' --outdir './'
```

Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json 'training_database.k31_ani0.50_config.json' --sample_file 'sample_well_80.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --out_filename 'k31_ani0.50_result.xlsx' --outdir './'
```


### Results
Decreasing ANI to 0.50 and using a ksize of 31, YACHT finds XYZ