import subprocess
from os.path import exists
import os
import pandas as pd
import shutil

def test_full_workflow():
    """
    Uses a random selection of genomes and a random metagenome sketch
    :return: None
    """
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))  # currently one level above ./tests
    test_dir = os.path.join(script_dir, 'tests')
    data_dir = os.path.join(test_dir, 'testdata')
    out_prefix = "20_genomes_trained"
    full_out_prefix = os.path.join(data_dir, out_prefix)
    abundance_file = os.path.join(data_dir, "result.xlsx")
    reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
    sample_sketches = os.path.join(data_dir, "sample.sig.zip")
    intermediate_dir = out_prefix + "_intermediate_files"
    # at this point, just checking one signature file
    expected_files = list(map(lambda x: os.path.join(data_dir, x), ["20_genomes_trained_config.json", "20_genomes_trained_processed_manifest.tsv"]))
    expected_files.extend(list(map(lambda x: os.path.join(data_dir, intermediate_dir, x), ["SOURMASH-MANIFEST.csv", "training_multisearch_result.csv",
                                              "training_sig_files.txt"])))
    # one of the signature files
    expected_files.extend(list(map(lambda x: os.path.join(data_dir, intermediate_dir, "signatures", x), ["04212e93c2172d4df49dc5d8c2973d8b.sig.gz"])))
    # remove the files if they exist
    for f in expected_files:
        if exists(f):
            os.remove(f)
    # Remove the intermediate folder
    shutil.rmtree(os.path.join(data_dir, intermediate_dir), ignore_errors=True)
    #  python ../make_training_data_from_sketches.py --ref_file testdata/20_genomes_sketches.zip --ksize 31 --prefix 20_genomes_trained --outdir testdata/
    cmd = f"python {os.path.join(script_dir, 'make_training_data_from_sketches.py')} --ref_file {reference_sketches}" \
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
    cmd = f"python {os.path.join(script_dir, 'run_YACHT.py')} --json {os.path.join(data_dir, '20_genomes_trained_config.json')} --sample_file {sample_sketches} --significance 0.99 --min_coverage 0.001 --out {os.path.join(data_dir,abundance_file)} --show_all"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output file exists
    assert exists(abundance_file)
    # check if CP032507.1 has correct results
    df = pd.read_excel(abundance_file)
    present_organism = "CP032507.1 Ectothiorhodospiraceae bacterium BW-2 chromosome, complete genome"
    # but not enough to claim presence
    assert str(df[df['organism_name'] == present_organism]["in_sample_est"].values[0]) == "False"
    # and we only observed 2 k-mers in the sample
    assert df[df['organism_name'] == present_organism]["num_matches"].values[0] == 2
    # and the threshold was 706
    assert df[df['organism_name'] == present_organism]["acceptance_threshold_with_coverage"].values[0] == 706


def test_incorrect_workflow1():
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    demo_dir = os.path.join(script_dir, "demo")
    cmd = f"python run_YACHT.py --json {demo_dir}/demo_ani_thresh_0.95_config.json --sample_file {demo_dir}/ref.sig.zip"
    res = subprocess.run(cmd, shell=True, check=False)
    # this should fail
    assert res.returncode == 1


def test_demo_workflow():
    cmd = "cd demo; sourmash sketch dna -f -p k=31,scaled=1000,abund -o sample.sig.zip query_data/query_data.fq"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; sourmash sketch fromfile ref_paths.csv -p dna,k=31,scaled=1000,abund -o ref.sig.zip --force-output-already-exists"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../make_training_data_from_sketches.py --force --ref_file ref.sig.zip --ksize 31 --num_threads 1 --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../run_YACHT.py --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads 1 --min_coverage_list 1 0.6 0.2 0.1 --out result.xlsx"
    _ = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../srcs/standardize_yacht_output.py --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0
    assert exists('./cami_result.cami')
