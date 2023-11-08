# Pathogen Detection Example
A patient with respiratory symptoms seeks to find out the pathogen that is causing them these symptoms.

Make sure all bacterial reads needed to create your reference dataset also known as a training dataset are available.
```bash
bash 1_before_starting.sh
```
```bash
bash 2_before_starting.sh
```

### Sketch your training dataset and sample to your preference.
Note: training and sample datasets are required to have the same ksize.

#### Using k=31
```bash
sourmash sketch fromfile reference_list.csv -p dna,k=31,scaled=1000,abund -o training_database.k31.sig.zip
```
### Make training data for k=31
```bash
python make_training_data_from_sketches.py --ref_file training_database.k31.sig.zip --ksize 31 --num_threads 32 --ani_thresh 0.95 --prefix 'training_database.k31.0.95' --outdir ./
```

### Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json '../training_database.k31.0.95_config.json' --sample_file 'lung_sample.k31.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --outdir './'
```

### Results
Using a ksize of 31, YACHT finds that M. pneumoniae is present in the lung sample.

## What if we decrease ksize to 15?
If we use small ksizes like 15, we would expect to not find that the patient is infected by M. pneumoniae. Let's set up the experiment

### Sketch Lung Sample using a k=15
```bash
sourmash sketch fromfile reference_list.csv -p dna,k=15,scaled=1000,abund -o training_database.k15.sig.zip
```
### Make training data for k=15
```bash
python make_training_data_from_sketches.py --ref_file training_database.k15.sig.zip --ksize 15 --num_threads 32 --ani_thresh 0.95 --prefix 'training_database.k15.0.95' --outdir ./
```
### Pathogen Detection using YACHT
Identify whether the patient has a infectin and what pathogen is causing the disease.
```bash
python ../../run_YACHT.py --json '../training_database.k15.0.95_config.json' --sample_file 'lung_sample.k15.sig.zip' --significance 0.99 --min_coverage 1 0.5 0.1 0.05 0.01 --outdir './'
```
### Results
Using a ksize of 15, YACHT finds/does not fine that M. pneumoniae