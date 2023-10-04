import subprocess
from os.path import exists
import os
import pandas as pd


def test_full_workflow():
    """
    Uses a random selection of genomes and a random metagenome sketch
    :return: None
    """
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))  # currently one level above ./tests
    test_dir = os.path.join(script_dir, 'tests')
    data_dir = os.path.join(test_dir, 'testdata')
    out_prefix = "integration_test"
    full_out_prefix = os.path.join(data_dir, out_prefix)
    abundance_file = full_out_prefix + "recovered_abundance.xlsx"
    reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
    sample_sketches = os.path.join(data_dir, "sample.sig")
    expected_files = map(lambda x: full_out_prefix + x, ["_hash_to_col_idx.csv", "_processed_org_idx.csv",
                                              "_ref_matrix_processed.npz", "_ref_matrix_unprocessed.npz",
                                              "_recover_abundance.csv", "_ksize_ani_thresh.json"])
    # remove the files if they exist
    for f in expected_files:
        if exists(f):
            os.remove(f)
    cmd = f"python {os.path.join(script_dir, 'make_training_data_from_sketches.py')} --ref_file {reference_sketches}" \
          f" --out_prefix {full_out_prefix} --ksize 31"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output files exist
    for f in expected_files:
        assert exists(f)
    # check that the files are big enough
    for f in expected_files:
        assert os.stat(f).st_size > 1000

    # then do the abundance estimation
    if exists(abundance_file):
        os.remove(abundance_file)
    cmd = f"python {os.path.join(script_dir, 'run_YACHT.py')} --json {os.path.join(script_dir, 'gtdb_ani_thresh_0.95_config.json')} --sample_file {sample_sketches} --significance 0.99 --min_coverage 1 --outdir {data_dir} --out_filename {abundance_file}"
    # print(cmd)
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output file exists
    assert exists(abundance_file)
    # check if CP032507.1 has correct abundance of 6
    df = pd.read_excel(abundance_file)
    present_organism = "CP032507.1 Ectothiorhodospiraceae bacterium BW-2 chromosome, complete genome"
    # but not enough to claim presence
    assert len(df[df['organism_name'] == present_organism]["in_sample_est"].values) == 0
    # and we only observed 2 k-mers in the sample
    assert len(df[df['organism_name'] == present_organism]["num_matches"].values) == 0
