# Produce Pathogen Detection Toy Dataset
# This script creates the toy datasets for pathogen detection use case example
# Make sure SRA tools is in your environment for these use case examples

## Bacterial reads
# As a naive type dataset would it matter where sample comes from?
### Download 20 reads from SRR25626360 which represents 20 reads from WGS of Haemophilus influenzae 
fastq-dump -X 20 -Z SRR25626360 > hinfluenzae_SRR25626360_20reads.fastq

### Download 20 reads from SRR24210460 which represents 20 reads from WGS of mycoplasma pneumoniae from library MDY 
fastq-dump -X 20 -Z SRR24210460 > mpneumoniae_SRR24210460_20reads.fastq

### Download 20 reads from SRR7217470 which represents 20 reads from WGS of Chlamydia pneumoniae 
fastq-dump -X 20 -Z SRR7217470 > xpneumoniae_SRR7217470_20reads.fastq

### Download 20 reads from SRR5962942 which represents 20 reads from WGS of Streptococcus pneumoniae 
fastq-dump -X 20 -Z SRR5962942 > spneumoniae_SRR5962942_20reads.fastq

### Download 20 reads from SRR26202532 which represents 20 reads from WGS of Bordetella pertussis
fastq-dump -X 20 -Z SRR26202532 > bpertussis_SRR26202532_20reads.fastq

## viral reads
### Download 20 reads from SRR26589836 sra which represents 20 reads from WGS of covid 
fastq-dump -X 20 -Z SRR26589836 > sars_cov_2_SRR26589836_20reads.fastq

# Download 100 reads from SRR2830253 which are reads of a healthy human lung microbiome
fastq-dump -X 100 -Z SRR2830253 > healthy_lung_SRR2830253_100reads.fastq

# download 100 reads from SRR13286708 which is the Lung bacterial microbiome of critical COVID-19
fastq-dump -X 100 -Z SRR13286708 > critical_COVID_patient_lungs_SRR13286708_100reads.fastq

# Create toy dataset that is tested on a healthy lung fastq reads

# Create toy dataset where fastq reads from either covid or h influenzae are included and use yacht to analyze

# Create a toy dataset that is a referrnce database of known organisms to trigger respiritory disease 