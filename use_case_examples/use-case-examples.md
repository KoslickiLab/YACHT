# Biological Application

Does YACHT provide qunatitative data?
How sensitive is Yacht in detecting viruses?
Could we use YACHT to identify the amount of host DNA?
Is sequencing depth an issue for YACHT?
Can YACHT be used to identify the amount of contamination that there is in a sample?

## Contamination Detection

Identifying contamination of samples is an important step for downstream analysis. 

Failure to detect sample contaminants can bias community diversity and strain sharing identification, which lead to false claims in research. Additionally, contaminant detection is a challenge for low-biomass data.

There are two types of contaminants that yacht can detect. Specifically, external and cross contaminations. External contaminations includes any microbial DNA that can come from a researcher's native microbiome, experimental kits, and surrounding areas. Cross contamination can come from DNA extractions, sequencing index switching, and sample bleeding. Further, cross-contamination can complicate contimanation detection.

Here, we explore how we can use YACHT for contaminant detection with the following tutorials:
* External Contamination Detection 
* Cross Contamination Detection 
* Contamination in Low-biomass Datasets 

### Contaminant Detection Examples

#### External Contamination Detection

Data:

Commands:

Interpretation:

#### Cross Contamination Detection

Data:

Commands:

Interpretation:

#### Contamination in Low-biomass Datasets

## Pathogenic Detection

Early and accurate pathogen detection is vital for rapid diagnostic and clinical intervention. Some microorganisms are hard to grow, cultivate (i.e. viruses), and/or take a long time to grow (i.e. mycobacteriums and mold). Additionally, the origin of the sample can decrease pathogen detection sensitivity such as purulent samples having human DNA noise. These challenges can lead to extended hospitalizations, readmissions, as well as increased mortality.

Research in pathogen detection have led to improvements of the experimental part such as differential lysis of human cells, etc. Further, efforts in improving computational tools and pipelines for pathogen detection have also been made such as real time pathogen detection (EPI2ME, SURPI-RT, etc)

Methods for pathogen detection include PCR, multiplex PCR, broad-range PCR, antigen detection, MALDI-TOF MS, and PNA-FISH. However, PCR based methods are limited to false negatives, knowledge priori, and/or viral detection. Further, antigen detection is limited to the time for antibodies to develop which can take 1-2 weeks. Finally, methods such as mass spectrometry and insitu hybridization have been advancements for bacterial and fungal identification but still require bacterial cltivation, do not provide quantitative results (in the case of MALDI-TOF MS) or are limited to tissue samples (in the case of PNA-FISH). 

Using clinical metagenomic data, we can detect any pathogenic infection (i.e. respiratory, bloodstream, central nervous system infections) without prior knowledge from millions of reads produced by shot gun metagenomic sequencing. 

### Pathogenic Detection Examples

