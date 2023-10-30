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

Lou YC, Hoff J, Olm MR, West-Roberts J, Diamond S, Firek BA, Morowitz MJ, Banfield JF. Using strain-resolved analysis to identify contamination in metagenomics data. Microbiome. 2023 Mar 2;11(1):36. doi: 10.1186/s40168-023-01477-2. PMID: 36864482; PMCID: PMC9979413.

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

Early and accurate pathogen detection is vital for rapid diagnostic, clinical intervention, pathogen discovery, pathogen surveillance, and outbreak investigations. Some microorganisms are hard to grow, cultivate (i.e. viruses), and/or take a long time to grow (i.e. mycobacteriums and mold). Additionally, the origin of the sample can decrease pathogen detection sensitivity such as purulent samples having human DNA noise. These challenges can lead to extended hospitalizations, readmissions, as well as increased mortality.

Research in pathogen detection have led to improvements of the experimental part such as differential lysis of human cells, etc. Further, efforts in improving computational tools and pipelines for pathogen detection have also been made such as real time pathogen detection (EPI2ME, SURPI-RT, etc)

Methods for pathogen detection include PCR, multiplex PCR, broad-range PCR, antigen detection, MALDI-TOF MS, and PNA-FISH. However, PCR based methods are limited to false negatives, knowledge priori, and/or viral detection. Further, antigen detection is limited to the time for antibodies to develop which can take 1-2 weeks. Finally, methods such as mass spectrometry and insitu hybridization have been advancements for bacterial and fungal identification but still require bacterial cultivation, do not provide quantitative results (in the case of MALDI-TOF MS) or are limited to tissue samples (in the case of PNA-FISH). 

Using clinical metagenomic data, we can detect any pathogenic infection (i.e. respiratory, bloodstream, central nervous system infections) without prior knowledge from millions of reads produced by shot gun metagenomic sequencing. However, pathogen detection can be difficult depending on the sample type being used such as purulent samples having a lot of human DNA nouse.

Computationally, pathogen detection is challenging due to alignment/classification algorithms becoming overwhelmed by large data, read sparsity (leading to de novo assembly difficulty), the lack of genome representation of novel pathogens, and the detection of pathogens that are highly divergent novel apthogens. Computationl subtraction is a popular approach to use for pathogen detection, tools such as PathSeq and SURPI use this approach. These approaches use alignments and although alignment based methods have been shown to slow down analyses, SURPI has shown to detect pathogens and generate results in approaite clinical timing from minutes to hours.

Gu W, Deng X, Lee M, Sucu YD, Arevalo S, Stryke D, Federman S, Gopez A, Reyes K, Zorn K, Sample H, Yu G, Ishpuniani G, Briggs B, Chow ED, Berger A, Wilson MR, Wang C, Hsu E, Miller S, DeRisi JL, Chiu CY. Rapid pathogen detection by metagenomic next-generation sequencing of infected body fluids. Nat Med. 2021 Jan;27(1):115-124. doi: 10.1038/s41591-020-1105-z. Epub 2020 Nov 9. PMID: 33169017; PMCID: PMC9020267.

Batool M, Galloway-Peña J. Clinical metagenomics-challenges and future prospects. Front Microbiol. 2023 Jun 28;14:1186424. doi: 10.3389/fmicb.2023.1186424. PMID: 37448579; PMCID: PMC10337830.

Naccache SN, Federman S, Veeraraghavan N, Zaharia M, Lee D, Samayoa E, Bouquet J, Greninger AL, Luk KC, Enge B, Wadford DA, Messenger SL, Genrich GL, Pellegrino K, Grard G, Leroy E, Schneider BS, Fair JN, Martínez MA, Isa P, Crump JA, DeRisi JL, Sittler T, Hackett J Jr, Miller S, Chiu CY. A cloud-compatible bioinformatics pipeline for ultrarapid pathogen identification from next-generation sequencing of clinical samples. Genome Res. 2014 Jul;24(7):1180-92. doi: 10.1101/gr.171934.113. Epub 2014 Jun 4. PMID: 24899342; PMCID: PMC4079973.

### Pathogenic Detection Examples

