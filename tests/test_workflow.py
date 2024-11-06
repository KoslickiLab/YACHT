import subprocess
from os.path import exists
import os
import pandas as pd
import shutil
import sys
cpath = os.path.dirname(os.path.realpath(__file__))
project_path = os.path.join(cpath,'..')

def test_full_workflow():
    """
    Uses a random selection of genomes and a random metagenome sketch
    :return: None
    """
    test_dir = os.path.join(project_path, 'tests')
    data_dir = os.path.join(test_dir, 'testdata')
    out_prefix = "20_genomes_trained"
    full_out_prefix = os.path.join(data_dir, out_prefix)
    abundance_file = os.path.join(data_dir, "result.xlsx")
    reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
    sample_sketches = os.path.join(data_dir, "sample.sig.zip")
    intermediate_dir = out_prefix + "_intermediate_files"
    # at this point, just checking one signature file
    expected_files = list(map(lambda x: os.path.join(data_dir, x), ["20_genomes_trained_config.json", "20_genomes_trained_processed_manifest.tsv"]))
    expected_files.extend(list(map(lambda x: os.path.join(data_dir, intermediate_dir, x), ["SOURMASH-MANIFEST.csv", "selected_result.tsv",
                                              "training_sig_files.tsv"])))
    # one of the signature files
    expected_files.extend(list(map(lambda x: os.path.join(data_dir, intermediate_dir, "signatures", x), ["04212e93c2172d4df49dc5d8c2973d8b.sig"])))
    # remove the files if they exist
    for f in expected_files:
        if exists(f):
            os.remove(f)
    # Remove the intermediate folder
    shutil.rmtree(os.path.join(data_dir, intermediate_dir), ignore_errors=True)
    #  python ../make_training_data_from_sketches.py --ref_file testdata/20_genomes_sketches.zip --ksize 31 --prefix 20_genomes_trained --outdir testdata/
    cmd = f"yacht train --ref_file {reference_sketches}" \
          f" --prefix {full_out_prefix} --ksize 31 --outdir {data_dir}"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output files exist
    for f in expected_files:
        assert exists(f)
    # check that the files are big enough
    for f in expected_files:
        assert os.stat(f).st_size > 291
    # then do the presence/absence estimation
    if exists(abundance_file):
        os.remove(abundance_file)
    # python ../run_YACHT.py --json testdata/20_genomes_trained_config.json --sample_file testdata/sample.sig.zip --out_file result.xlsx
    cmd = f"yacht run --json {os.path.join(data_dir, '20_genomes_trained_config.json')} --sample_file {sample_sketches} --significance 0.99 --min_coverage 0.001 --out {os.path.join(data_dir,abundance_file)} --show_all"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output file exists
    assert exists(abundance_file)
    # check if CP032507.1 has correct results
    df = pd.read_excel(abundance_file)
    present_organism = "CP032507.1 Ectothiorhodospiraceae bacterium BW-2 chromosome, complete genome"
    # but not enough to claim presence
    assert str(df[df['organism_name'] == present_organism]["in_sample_est"].values[0]) == "True"
    # and we only observed 2 k-mers in the sample
    assert df[df['organism_name'] == present_organism]["num_matches"].values[0] == 2
    # and the threshold was 706
    assert df[df['organism_name'] == present_organism]["acceptance_threshold_with_coverage"].values[0] == 0


def test_incorrect_workflow1():
    demo_dir = os.path.join(project_path, "demo")
    cmd = f"yacht run --json {demo_dir}/demo_ani_thresh_0.95_config.json --sample_file {demo_dir}/ref.sig.zip"
    res = subprocess.run(cmd, shell=True, check=False)
    # this should fail
    assert res.returncode == 1


def test_demo_workflow():
    cmd = f"cd {project_path}/demo; yacht sketch sample --infile ./query_data/query_data.fq --kmer 31 --scaled 1000 --outfile sample.sig.zip"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = f"cd {project_path}/demo; yacht sketch ref --infile ./ref_genomes --kmer 31 --scaled 1000 --outfile ref.sig.zip"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = f"cd {project_path}/demo; yacht train --force --ref_file ref.sig.zip --ksize 31 --num_threads 1 --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = f"cd {project_path}/demo; yacht run --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads 1 --min_coverage_list 1 0.6 0.2 0.1 --out result.xlsx"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = f"cd {project_path}/demo; yacht convert --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./"
    _ = subprocess.run(cmd, shell=True, check=True)


