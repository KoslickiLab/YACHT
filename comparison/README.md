# Benchmarking Experiments

To evaluate the performance of YACHT, we utilize the robust [public datasets](https://data.cami-challenge.org/participate) from the [Critical Assessment of Metagenome Interpretation II (CAMI II)](https://www.nature.com/articles/s41592-022-01431-4) for the tasks of taxonomic profiling and clinical pathogen detection, which contains the rhizosphere data, marine data, strain madness data, and clinical pathogen data. These datasets include the rhizosphere data, marine data, strain madness data, and clinical pathogen data, providing high-quality benchmarks for evaluating metagenomic analysis tools like YACHT. We also leverage the CAMI-official profiling assessment tool [OPAL](https://github.com/CAMI-challenge/OPAL/tree/master) to compare YACHT's performance against the state-of-the-art (SOTA) tools (e.g., Bracken, Metaglin, mOTUs, MetaPhlAn, CCMetagen, NBC++, MetaPhyler, LSHVec) suggested by CAMI II. 

## How to reproduce the evaluation results
We provide two bash scripts under `scripts` folder to reproduce our evaluation results. To run these scripts, please make sure your have installed the Conda environnment suggested in the repository home page. After that, run the following commands:
1. Download CAMI2 datasets 
```bash
# download_cami2_data.sh <output_dir> <cpu_num>
bash download_cami2_data.sh ~/YACHT/comparison 50
```
2. Run benchmarking experiments
```bash
run_experiments.sh <yacht_repo_loc> <output_dir> <cpu_num>
bash download_cami2_data.sh ~/YACHT ~/YACHT/comparison 20
```


