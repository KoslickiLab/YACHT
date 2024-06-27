# Create the example sample data for a patient with respiratory symptoms seeks to find out the pathogen that is causing them these symptoms.

# Before moving on. Make sure reads needed to create sample dataset are available. Please reference create_reference_database.md

# Sketch sample to your preference. Note: training and sample datasets are required to have the same ksize.

## Using k=31
nohup sourmash sketch fromfile lung_list.csv -p dna,k=31,scaled=1000,abund -o lung_sample.k31.sig.zip > k31_sample.log 2>&1 &

## Using k=15
nohup sourmash sketch fromfile lung_list.csv -p dna,k=15,scaled=1000,abund -o lung_sample.k15.sig.zip > k15_sample.log 2>&1 &
