# Pathogen Detection Example

### Pathogen detection create a lung sample with a spiked pathogen
Make sure SRA tools is in your environment for these use case examples

##### Use nohup to download reads from SRR2830253, reads of a healthy human lung microbiome
```bash
nohup fastq-dump --fasta 60 SRR2830253 2>&1 &
```

```bash
mv SRR2830253.fasta healthy_lung_SRR2830253.fa
```
```bash
nohup sourmash sketch dna -f healthy_lung_SRR2830253.fasta -p k=31,scaled=1000,abund -o healthy_lung_SRR2830253.sig.zip 2>&1 &
```

##### Use nohup to download reads from SRR24210460, reads WGS of mycoplasma pneumoniae from library MDY 
```bash
nohup fastq-dump --fasta 60 SRR24210460 2>&1 &
```

```bash
mv SRR24210460.fasta mpneumoniae_SRR24210460.fasta
```
```bash
nohup sourmash sketch dna -f mpneumoniae_SRR24210460.fasta -p k=31,scaled=1000,abund -o mpneumoniae_SRR24210460.sig.zip 2>&1 &
```

##### Merge signatures of mycoplasma pneumoniae and healthy lung to creat a toy dataset that is a human lung sample with mycoplasma pneumoniae
```bash
nohup sourmash signature merge healthy_lung_SRR2830253.sig.zip mpneumoniae_SRR24210460.sig.zip -o lung_sample.sig.zip > log 2>&1 &
```
