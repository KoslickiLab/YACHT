---
title: 'YACHT: an ANI-based statistical test to detect microbial presence/absence in a metagenomic sample'
tags:
  - Python
  - c++
  - metagenomics
  - microbial
authors:
  - name: Chunyu Ma
    orcid: 0000-0001-9731-5153
    equal-contrib: true
    affiliation: 1
  - name: Adam Park
    orcid: 0009-0003-4104-2960
    equal-contrib: true
    affiliation: 1
  - name: Shaopeng Liu 
    orcid: 0000-0003-3112-4068
    equal-contrib: true
    affiliation: 2
  - name: Omar Hesham Rady
    orcid: 0009-0005-1819-7643
    equal-contrib: true
    affiliation: 1
  - name: Mahmudur Rahman Hera
    orcid: 0000-0002-5992-9012
    equal-contrib: true
    affiliation: 1
  - name: Judith S. Rodriguez
    orcid: 0000-0002-5109-3054
    equal-contrib: true
    affiliation: 2
  - name: Stephanie J. Won
    orcid: 0000-0001-7288-5395
    equal-contrib: true
    affiliation: 3
  - name: Maksym Lupei
    orcid: 0000-0003-3440-3919
    equal-contrib: true
    affiliation: 1
  - name: David Koslicki
    corresponding: true 
    affiliation: "1, 2, 3"
affiliations:
  - name: School of Electrical Engineering and Computer Science, Pennsylvania State University, USA
    index: 1
    ror: 04p491231
  - name: Huck Institutes of the Life Sciences, Pennsylvania State University, USA
    index: 2
    ror: 04p491231
  - name: Department of Biology, Pennsylvania State University, USA
    index: 3
    ror: 04p491231
date: 10 February 2025
bibliography: paper.bib
---

# Summary

In metagenomics, identifying genomes present in a sample is an important initial task, but is complicated by taxonomic profiling tools lacking uncertainty quantification and using incomplete reference databases missing exact genome matches. YACHT (**Y**es/No **A**nswers to **C**ommunity membership via **H**ypothesis **T**esting) [@koslicki2024yacht] is a command-line tool for taxonomic profiling that uses hypothesis testing to confidently determine genome presence/absence in a metagenomic sample. YACHT assists in discovering rare microbiomes by identifying low-abundant species missed in other taxonomic profiling approaches while also controlling the false negative rate. Its statistical model overcomes challenges in sequencing coverage and incomplete genomes, making it ideal for diverse metagenomic applications, including functional profiling, metatranscriptomics, and clinical microbiome analysis.

YACHT presents a robust, $k$-mer sketching-based statistical framework for accurately detecting genetic similarity between the reference database and the metagenomic sample by incorporating evolutionary sequence divergence through the average nucleotide identity (ANI) and sequencing coverage to enable efficient detection of sampled genomes. The workflow for YACHT includes the following commands. To begin, `yacht sketch` creates reduced representation "sketches" of the reference and sample datasets enabling swift comparisons. Then, `yacht train` is used to find a representative of closely related reference genomes using ANI. Lastly, `yacht run` uses the YACHT algorithm to perform hypothesis testing and identify the presence or absence of species. YACHT is developed with C++ and Python and depends on `sourmash` [@irber2024sourmash], a program for extracting and managing $k$-mers.

# Statement of need

Accurately identifying and characterizing microbial communities with low relative abundance is a significant challenge in metagenomics. The current profiling-based practice involves setting arbitrary filter thresholds or discarding low-abundance data without robust justification, which can compromise profiling accuracy and lead to misinterpretations [@Schloss; @jia2022sequencing]. Even with such filtering, the results remain inherently arbitrary because they are influenced by biological complexities such as sequencing errors and evolutionary processes. The lack of a systematic approach to establishing credibility in these results diminishes researchers' confidence in biologically informed methods for identifying rare microorganisms, thereby undermining metagenomic studies. Moreover, these difficulties are exacerbated by the incompleteness of reference databases and the variability in sequencing coverage depth, underscoring the need for statistically credible approaches.

Metagenomic methods rely on existing genome references to detect and classify microbial organisms. However, these reference databases are often incomplete, and conventional metrics may not always align with traditional taxonomic frameworks that account for genomic changes. Consequently, microbes that carry mutations or have diverged evolutionarily can remain undetected, causing inaccuracies in microbial community profiling and misinterpretation of data [@kunin2008bioinformatician; @schlaberg2017validation; @loeffler2020improving; @marcelino2020ccmetagen]. Hence, analytical frameworks need to incorporate genome similarity metrics to capture the full breadth of microbial diversity and to provide accurate, interpretable microbiome dynamics. However, incomplete databases alone do not account for all metagenomic challenges; sequence coverage depth also contributes to the resolution and reliability of microbial detection and characterization.

Sequence coverage depth, defined as the portion of a microbe’s genome detected in a sample, is crucial for detecting low-abundance microbes. However, sequencing processes often fail to achieve complete coverage of all genomes in a sample due to limited sequencing depth. As a result, rare or low-abundance taxa may exhibit low sequence coverage, leading to their misinterpretation as noise rather than genuine observations [@mande2012classification; @shakya2013comparative; @sczyrba2017critical; @meyer2022critical]. Furthermore, the lack of guidelines for establishing a biologically meaningful coverage depth threshold introduces subjectivity and inconsistency in the metagenomic analyses. Therefore, implementing dynamic coverage depth thresholds tailored to varying abundance levels is essential for delivering accurate metagenomic studies. Yet, even if we address coverage depth and incomplete genome reference problems, ensuring proper control over statistical errors remains another major challenge.

Existing metagenomic methods lack the statistical rigor to control false positives and false negatives effectively. High false positive rates misrepresent microbial composition and lead to biased conclusions, undermining research reliability. Conversely, false negative rates cause researchers to overlook important taxa, especially those in low abundance that often carry significant biological importance [@jousset2017less]. Incomplete reference databases, sequencing errors, and evolutionary divergence between reference and sample genomes further complicate these challenges. Therefore, maintaining appropriate control over these statistical error rates is critical to ensure more confident, reliable biological inferences and minimize the risk of misinterpretation. While limitations in reference database, sequence coverage depth and balance of statistical error pose significant challenges, the complexity of metagenomic analysis demands a multifaceted approach to capture microbial profiling accurately.

To address these challenges, YACHT offers a statistical framework that can accurately determine the presence or absence of microbial genome in a sample through hypothesis testing. The algorithm’s mathematical model accounts for evolutionary sequence divergence and incomplete sequencing depth by utilizing genome similarity and minimum sequencing depth parameters. It employs the FrachMinHash sketching technique, an alignment-free $k$-mer approach, facilitating fast and accurate genome detection that can efficiently process large datasets. YACHT ensures precise detection of low abundance taxa with a user-defined false negative rate, minimizing the risk of misinterpretation of the result. Our approach can be used for other metagenomic applications such as functional profiling, metatranscriptomic studies [@marcelino2019metatranscriptomics], metabolic potential analyses [@ward2018metapoap; @pereira2024metatranscriptomics], and the characterization of low abundant clinical metagenomic samples such as skin [@godlewska2020metagenomic]. YACHT enhances metagenomic analysis by offering reduced reliance on arbitrary thresholds, improving the interpretability of the result without compromising biological relevance, and allowing researchers to differentiate between genuine artifacts from “noise” with statistical confidence.

![the workflow of the YACHT algorithm for deciding inclusion of microbial genome in samples. The process involves four main steps: Generate k-mer sketches for both reference genomes and samples using sourmash; Preprocess the reference genomes using the yacht train command to merge redundant genomes based on ANI thresholds, retaining the most representative genome; Run the YACHT algorithm to perform probabilistic hypothesis testing with user-defined parameters.](workflow.png)

# Workflow

## 1. Create sketches

This step uses the command `yacht sketch`, which generates sketches for both the samples and the reference genomes. Users must utilize sourmash to extract sketches from a reference database of microbial genomes. [sourmash Databases](https://sourmash.readthedocs.io/en/latest/databases.html) provide a variety of pre-formed databases of such sketches, or users can create a custom database using the `sourmash sketch` command on FASTA/FASTQ files of reference genomes (see the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/)). Other available databases include the [GTDB genomic representatives database](https://gtdb.ecogenomic.org/downloads). The sketches for samples must be generated using the same $k$-mer size and scale factor as those used for the reference database. The scale factor acts as an indicator of data compression, with smaller values being more appropriate for smaller datasets.

## 2. Preprocess the reference genomes

This step uses the command `yacht train`, which identifies and merges genomes that are roughly identical based on Average Nucleotide Identity (ANI). It extracts sketches from the reference database and converts them into a binary format specific to YACHT. For any two organisms with ANI greater than a user-specified threshold, one will be removed from the reference set, as such organisms are considered too similar to be "distinguishable."

For instance, if the threshold `ani_thresh` is set to 0.95, two organisms with ANI greater than 0.95 will be deemed indistinguishable. In such cases, the organism with the largest number of unique $k$-mers will be retained in the reference set. If there is a tie in the number of unique $k$-mers, one organism will be randomly selected.

## 3. Run the YACHT algorithm

This step uses the command `yacht run`, which executes the hypothesis test implementing the YACHT algorithm. Two parameters are of particular importance for controlling statistical power:

- `significance`: This parameter is akin to the confidence level, controlling how certain we want to be that the organism is present. A higher value generally leads to more false negatives, while a lower value leads to more false positives.
- `min_coverage_list`: This parameter specifies a list of min coverage values from 0 to 1 that indicate the percentage of distinct $k$-mers that must be sequenced and present in the sample to qualify an organism as being present. Setting this to 1 is a conservative approach, but if the sample has low coverage, a lower value may be more appropriate. A higher value will increase false negatives, while a lower value increases false positives. In practice, we find this parameter to be the most important for controlling the trade-off between precision and recall.

The outputs provide probabilistic decisions regarding the presence and absence of organisms, described as the following representative output columns. The column `num_matches` indicates the number of $k$-mers found in both the organism and the sample. The column `acceptance_threshold_*` specify the number of $k$-mers required to be found in the organism to consider it present at the given ANI threshold as indicated by either TRUE or FALSE in the `Presence_*` column. For example, the presence for organisms *Sediminispirochaeta* and *Natronobacterium* are reported as TRUE because the total $k$-mers matched to the reference is equal or greater than the required total matches indicated by the `acceptance_threshold_*` column. In contrast, the organism *Echinicola* does not have the required $k$-mer matches to be considered as present, therefore its presence is reported as FALSE. Finally, the columns `alt_confidence_mut_rate_*` denote the mutation rate (1-ANI) necessary for the false positive rate to match the false negative rate of 1-`significance`. I.e. if one minus this value is the true ANI of genome in the sample to the reference genome, then the probability of a false positive is the same as the false negative rate (based on significance). So the closer this value is to the given ANI threshold, the greater the statistical power.


| Organism             | Presence | num_matches | acceptance_threshold | alt_confidence_mut_rate |
|----------------------|----------|----------------|-------------------------|----------------------------|
| Sediminispirochaeta  | TRUE     | 2572           | 895                     | 0.053008659                |
| Natronobacterium     | TRUE     | 700            | 638                     | 0.053534755                |
| Echinicola           | FALSE    | 244            | 978                     | 0.052885411                |

**Table:** Example results with an ANI threshold of 0.95 and a minimum coverage set to 1. The column `num_matches` reports the number of $k$-mers shared between the organism and the sample. The `acceptance_threshold` indicates how many $k$-mers must be observed in the organism to consider it present at the specified ANI threshold. Finally, `alt_confidence_mut_rate` reflects the mutation rate (1 - ANI) at which the false positive rate would equal the false negative rate, defined as (1 - the significance level).

## 4. Convert the results to other popular output formats
This step uses the command `yacht convert`, which converts the results into other popular output formats, including CAMI, BIOM, and GraphPhlAn. For the conversion, a TSV file containing two columns, genome ID and its corresponding taxid, must be provided.


## Remarks on performance
The YACHT algorithm is lightweight as it operates on reduced datasets of $k$-mers, i.e., sketches. However, the second step of the workflow—preprocessing the reference genomes—can require significant time to identify indistinguishable sets of references when implemented with sourmash (keeping in mind that this step needs only be performed once on a given reference database and the results can be used to analyze multiple samples). To address this, we have developed a C++ implementation to reduce the computational overhead and ensure that the second step does not become a bottleneck.


The `yacht train` command performs pairwise computations on sequences, merging "indistinguishable" references based on ANI to generate a compact reference set. In our test using the GTDB representative genomes (r214), which includes 85,205 species-level genomes, `yacht train` required approximately 12 minutes and 52 GB of RAM to preprocess the data on an Ubuntu 22.04.5 system with 64 threads. If only a subset of genomes is needed, the `sourmash sig` command can be used to extract the signatures of interest for comparison.

# Use case examples

We present several use case examples demonstrating the application of YACHT for identifying taxonomy in microbiome studies. Specifically, we focus on three applicable scenarios that include low abundance samples, metagenome-assembled genome (MAG) fishing, and synthetic metagenomes. With these use case examples, we explore YACHT's flexible parameters, such as $k$-size, ANI, and minimum coverage, demonstrating how modifying these settings can influence and yield valuable insights into the taxonomy of these samples.

## Case I: Low abundance samples

Metagenomic samples, such as those obtained from skin, can yield low microbial DNA concentrations, which can lead to biased conclusions. YACHT provides users with the flexibility to modify key parameters, such as the ANI threshold and $k$-size, to enhance the detection of difficult-to-identify organisms. Thus, adjusting these parameters can improve reliability in species identification. For example, increasing the ANI threshold restricts species merging in preprocessing, while increasing the $k$-size makes YACHT more stringent when comparing samples to the trained reference. In this use case example, we demonstrate how tuning these parameters influences species reporting.

We applied YACHT to a randomly selected metagenomic sample, SRR32008482 (Bioproject PRJNA1209946)—a human facial skin sample sequenced from whole shotgun sequencing. For the following discussion, we focus our attention to a moderate minimal coverage of 0.05. By varying $k$-sizes (21, 31, and 51) and ANI thresholds (0.90, 0.95, and 0.995), we analyzed how these adjustments impact genome reporting. First, at an ANI threshold of 0.95, we observed that increasing the $k$-size from 21 to 51 resulted in more unique species being reported—rising from 19 species at $k$-size 21 to 115 species at $k$-size 51. Importantly, larger $k$-sizes experience a trade-off in that more exclusive $k$-mers can be identified corresponding to reference genomes, but higher sensitivity to sequencing error and mutations, resulting in few matches to a sample (To ascertain the “correct” $k$-mer size to use, one can fix all other parameters, start with $k=31$, vary the $k$-mer size up and down, and choose the size where the total number of predictions is minimized). Consequently, a higher $k$-size would require the sample to have more $k$-mers matching to the genome to confidently report as present. Too high of a $k$-mer size will start introducing false positives. Second, adjusting the ANI threshold while using a $k$-size of 31, we found that an ANI of 0.995 enables strain-level resolution by distinguishing highly similar but distinct strains. Conversely, at lower ANI thresholds, species differentiation becomes less precise. When using an ANI of 0.90, YACHT may report species that are broadly similar, thus careful interpretation is needed. To mitigate this, we recommend increasing the $k$-size when using a lower ANI threshold to improve specificity.

## Case II: Metagenomic-assembled genome fishing

In the case of metagenomic-assembled genome (MAG) fishing, we can employ YACHT to search or “fish” for a MAG of interest. In the previous case, we used sample SRR32008482 to demonstrate how modifying parameters can influence species reporting. Taking a closer look at that analysis, we observe the detection of *Bacteroidales bacterium* (Accession: GCA_017506175) in sample SRR32008482 when using a $k$-size of 51 and an ANI threshold of 0.95 at a minimum coverage of 0.5. To elaborate, *B. bacterium* has 2,464 exclusive $k$-mers, and given our parameters, a minimum of 158 matched $k$-mers was necessary to confidently identify it as present. For sample SRR32008482, there were 1,232 matched $k$-mers, clearly surpassing the threshold requirement to report as present. Expanding on the usage of YACHT, we can perform MAG Fishing—employing a single MAG as the training reference database rather than an extensive database such as GTDB. To illustrate MAG fishing using YACHT, we selected a second metagenomic sample, SRR32008483 (also from Bioproject PRJNA1209946), to assess the presence of *B. bacterium* using similar parameters (i.e., $k$-size=51, ANI=0.95, and scale factor=10). When employing a single MAG in YACHT, signatures are generated for every entry found in the FASTA/Q files of the MAG of interest. Therefore, using *B. bacterium* as a MAG reference database, a signature is produced for every scaffold entry, resulting in the detection of *B. bacterium* scaffolds in SRR32008483 at low minimum coverages of 0.1, 0.05, and 0.01. In contrast, for the previous sample SRR32008482, *B. bacterium* scaffolds are present across a broader range of minimum coverages (1.0, 0.5, 0.1, 0.05, and 0.01). The variation in matched $k$-mers across these different coverage thresholds between samples SRR32008482 and SRR32008483 highlights the critical influence of sequencing coverage. Specifically, SRR32008483 has a notably lower sequencing coverage (6.6 Mb) compared to SRR32008482 (44 Mb), reducing the number of exclusive $k$-mers available to match the MAG of interest. Additionally, comparing scaffolds shared between these samples at varying coverage thresholds, we observe that sample SRR32008482 contains 41 unique scaffolds of *B. bacterium* at a minimum coverage of 0.5 that are not present in sample SRR32008483 at a lower minimum coverage of 0.1. Overall, these observations demonstrate the impact of sequencing depth on MAG fishing.

## Case III: Synthetic metagenomes

Mock or synthetic microbe communities—also known as syncoms—serve as reference standards that are used in metagenomics to benchmark and optimize sequencing-based methods [@hardwick2018synthetic; @singer2016next], as well as evaluating the taxonomic and functional diversity of native microbial communities [@van2023synthetic]. YACHT can assess the construction of these synthetic communities to verify the designed microbes are present in the mock community for method evaluation and experimental validation. This can be helpful for identifying organisms missing from (or unexpectedly present in) a synthetic community [@tian2025designed].

Previously, we explored how $k$-size and minimum coverage can influence results made by YACHT. Here, we focus on ANI and false positives. Note carefully that we are not comparing at the species level, rather requiring same genome/strain matching for presence and absence. For this use case example, we use an SRA sample (SRR6394747) from the Bioproject PRJNA422917, a mock community called MBARC-26 composed of 26 microbes spanning bacteria and archaea species generated from whole shot-gun sequencing of long reads [@singer2016next]. Given the prior knowledge of the genomes expected in this syncom, we exclude five microbes, given that they are not found in the reference training database, and evaluate the ability of YACHT to identify the presence of 21 genomes using a $k$-size of 31 and minimum coverage of 0.05 across various ANIs of 0.80, 0.95, and 0.9995.

At an ANI of 0.95, YACHT identifies 19 expected genomes, with only 8 of those not being reported as present, potentially due to lower sequencing coverages resulting in a reduction of exclusive $k$-mers. When ANI is restricted to 0.9995, YACHT identifies 16 expected genomes barely making a difference when using an ANI of 0.95. On the other hand, if we loosen the constraint of ANI to 0.80, YACHT, as expected, produces many false positives, identifying only 3 of the expected genomes. This is due to many reference organisms (distinct from the genomes in the sample) sharing an ANI of 0.80 with those in the sample. As we lower the ANI from 0.80, an increasing number of false positives are reported as present. If we dramatically reduce the ANI, we will observe a substantial uptick in species reported as present, indicating an overly permissive ANI threshold will inflate false positives. As previously noted in the “Low abundance samples” case, species differentiation becomes less reliable, warranting precaution when using low ANIs such as 0.80, as 80% ANI is a very lenient definition of “same genome.” We also examined the effect of increasing the minimum coverage threshold to 0.25 while using an ANI of 0.95. In this case, YACHT reports fewer genomes as present, reducing false positives but also filtering out some true positives—highlighting the tradeoff between sensitivity and specificity at higher coverage thresholds.

Although we focused on the 21 expected genomes from the synthetic community, this is not a restriction imposed by YACHT, which considers all organisms in the reference database. When we broaden our evaluation to include all genomes reported by YACHT, the patterns remain consistent. At ANI 0.95 and coverage 0.05, YACHT identifies 23 genomes, 22 of which match the sample. At a stricter ANI of 0.9995 with the same coverage, 20 genomes are reported, but only 13 match. At a lenient ANI of 0.80, YACHT detects 74 genomes, with just 3 corresponding to the sample—suggesting the mock community likely comprises a mix of exact strain matches and closely related isolates. Some genomes in the community may not be perfectly represented in the NCBI reference set, but instead belong to the same species as the available reference genomes.

# Acknowledgements
We thank the contributors and collaborators who supported the development of YACHT. This work was supported in part by the National Institutes of Health (NIH) under grant number 5R01GM146462-03.

# References