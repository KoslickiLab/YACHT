import pandas as pd
import zipfile
import glob
import sys
import os
import yacht
import sourmash
from multiprocessing import Pool
import multiprocessing
from yacht.hypothesis_recovery_src import get_exclusive_hashes
from yacht.hypothesis_recovery_src import get_organisms_with_nonzero_overlap
from yacht.hypothesis_recovery_src import hypothesis_recovery
from yacht.utils import decompress_all_sig_files
from yacht.utils import load_signature_with_ksize
from yacht.utils import load_one_sig
"""
A script that reproduces the results of hypothesis_recovery_src so that the sylph coverage model can be incorporated
"""
manifestFilePath="/Users/randolph.raborn/Desktop/YACHT/demo/query_data/gtdb-rs214-reps.k31_0.95_pretrained/gtdb-rs214-reps.k31_0.95_processed_manifest.tsv"
sampleSigPath="/Users/randolph.raborn/Desktop/YACHT/demo/new_demo_files/SRR6940089_sample_out.sig.zip"
genTempDir="/Users/randolph.raborn/Desktop/YACHT/demo/query_data/gtdb_ani_thresh_0.95_intermediate_files/"
sampleTempDir="/Users/randolph.raborn/Desktop/YACHT/demo/new_demo_files/genome_reference.sig_temp"
multisearch_result_file = "/Users/randolph.raborn/Desktop/YACHT/demo/query_data/gtdb-rs214-reps.k31_0.95_pretrained/gtdb-rs214-reps.k31_0.95_intermediate_files/training_multisearch_result.csv"
multisearchResultPath= "/Users/randolph.raborn/Desktop/YACHT/demo/new_demo_files/genome_reference.sig_temp/sample_multisearch_result.csv"
ksize=31
min_cov=0.10 #setting to 0.10

sample_sig = load_signature_with_ksize(sampleSigPath, ksize)
print(sample_sig)
sample_info_set = (sampleSigPath, sample_sig)

#importing the manifest: 
manifest = pd.read_csv(manifestFilePath, sep="\t", header=0)
#print(manifest.head())

#following the procedure at the end of yacht's `get_organisms_with_nonzero_overlap`
#multisearch_result = pd.read_csv(
#            multisearch_result_file,
#            sep=",",
#            header=0,
#        )

#print(multisearch_result)
#multisearch_result_new = multisearch_result.drop_duplicates().reset_index(drop=True) 
#print(multisearch_result_new)
#multisearch_result_names = multisearch_result["match_name"].to_list() #this is what is actually returned by get_organisms_with_nonzero_overlap
#print(set(multisearch_result_names)) 

multiprocessing.set_start_method('fork')
#multisearch_result_new2 = get_organisms_with_nonzero_overlap(manifest, sampleSigPath, 1000, 31, 2, genTempDir, sampleTempDir) #comment this out to import directly from the multisearch result file
#print(multisearch_result_new2)
#multisearch_result_new2_file = pd.read_csv(multisearchResultPath, sep=",", header=0)
#multisearch_result_new2_file2 = multisearch_result_new2_file.drop_duplicates().reset_index(drop=True)

#print(f"Manifest")
#print(manifest.head())
#print(f"Multisearch_object")
#print(multisearch_result_new2_file2.head())

#multisearch_result_names2 = multisearch_result_new2_file["match_name"].to_list()
print(f"Made it here")
#print(multisearch_result_names2)
#exclusive_hashes_out = get_exclusive_hashes(manifest, multisearch_result_names2, sample_sig, ksize, genTempDir)
#print(exclusive_hashes_out)

#creating a minimum coverage list
min_cov_list = [min_cov] * len(manifest)

hyp_rec_out = hypothesis_recovery(manifest=manifest, sample_info_set=sample_info_set, path_to_genome_temp_dir=genTempDir, min_coverage_list=min_cov_list, scale=1000, ksize=ksize, ani_thresh=0.85, num_threads=16)

print(hyp_rec_out)