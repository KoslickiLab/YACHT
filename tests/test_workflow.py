import subprocess
from os.path import exists
import tempfile
import numpy as np
import os
from scipy import sparse
import pandas as pd


def test_small_edge_lengths():
    """
    Uses a random selection of genomes and a random metagenome sketch
    :return: None
    """
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))  # currently one level above ./tests
    data_dir = "testdata"
    out_prefix = "unittest_"
    full_out_prefix = os.path.join(data_dir, out_prefix)
    abundance_file = full_out_prefix + "recovered_abundance.csv"
    reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
    sample_sketches = os.path.join(data_dir, "sample.sig")
    expected_files = map(lambda x: full_out_prefix + x, ["hash_to_col_idx.csv", "processed_org_idx.csv",
                                              "ref_matrix_processed.npz", "ref_matrix_unprocessed.npz"])
    # remove the files if they exist
    for f in expected_files:
        if exists(f):
            os.remove(f)
    cmd = f"python {os.path.join(script_dir, 'ref_matrix.py')} --ref_file {reference_sketches} --out_prefix" \
          f" {full_out_prefix} --ksize 31"
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
    cmd = f"python {os.path.join(script_dir, 'recover_abundance.py')} --ref_file {full_out_prefix}ref_matrix_processed.npz --sample_file " \
          f"{sample_sketches} --w 0.01 --outfile {abundance_file} --ksize 31"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output file exists
    assert exists(abundance_file)
    # check if all the abundances are zero
    df = pd.read_csv(abundance_file, sep=",", header=0)
    assert np.allclose(df["estimated abundance"].values, np.zeros_like(df["estimated abundance"].values), atol=1e-2)
    # then run it again with a different w
    if exists(abundance_file):
        os.remove(abundance_file)
    cmd = f"python {os.path.join(script_dir, 'recover_abundance.py')} --ref_file {full_out_prefix}ref_matrix_processed.npz --sample_file " \
          f"{sample_sketches} --w 0.0001 --outfile {abundance_file} --ksize 31"
    #print(cmd)
    res = subprocess.run(cmd, shell=True, check=True)
    # check that no errors were raised
    assert res.returncode == 0
    # check that the output file exists
    assert exists(abundance_file)
    # check if all the abundances are zero
    df = pd.read_csv(abundance_file, sep=",", header=0)
    assert df[df['organism name'] == "CP032507.1 Ectothiorhodospiraceae bacterium BW-2 chromosome, complete genome"][
               "estimated abundance"].values[0] == 6.0
