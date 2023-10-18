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
    out_prefix = "20_genomes_trained"
    full_out_prefix = os.path.join(data_dir, out_prefix)
    abundance_file = os.path.join(data_dir, "result.xlsx")
    reference_sketches = os.path.join(data_dir, "20_genomes_sketches.zip")
    sample_sketches = os.path.join(data_dir, "sample.sig.zip")
    expected_files = map(lambda x: full_out_prefix + x, ["_hash_to_col_idx.csv", "_processed_org_idx.csv",
                                              "_ref_matrix_processed.npz", "_ref_matrix_unprocessed.npz",
                                              "_recover_abundance.csv", "_ksize_ani_thresh.json"])
    # In test data
    # 20_genomes_trained_config.json
    # 20_genomes_trained_processed_manifest.tsv
    # 20_genomes_trained_rep_to_corr_orgas_mapping.tsv
    # /testdata/20_genomes_trained_intermediate_files$ tree .
    # .
    # ├── signatures
    # │   ├── 04212e93c2172d4df49dc5d8c2973d8b.sig.gz
    # │   ├── 04f2b0e94f2d1f1f5b8355114b70274e.sig.gz
    # │   ├── 0661ecab88c3d65d0f10e599a5ba1654.sig.gz
    # │   ├── 06ebe48d527882bfa9505aba8e31ae23.sig.gz
    # │   ├── 11fe9a00287c7ad086ebbc463724cf10.sig.gz
    # │   ├── 16c6c1d37259d83088ab3a4f5b691631.sig.gz
    # │   ├── 188d55801a78d4773cf6c0b46bca96ba.sig.gz
    # │   ├── 1a121dca600c6504e88252e81004f0cf.sig.gz
    # │   ├── 39ea7fd48ee7003587c9c763946d5d6e.sig.gz
    # │   ├── 45f2675c1dca4ef1a24a05f5b268adbb.sig.gz
    # │   ├── 7b312fffa3fb35440ba40203ba826c05.sig.gz
    # │   ├── 8fb9b1a69838a58cc4f31c1e42a5f189.sig.gz
    # │   ├── 92fb1b3e4baa6c474aff3efb84957687.sig.gz
    # │   ├── 96cb85214535b0f9723a6abc17097821.sig.gz
    # │   ├── a136145bee08846ed94c0406df3da2d4.sig.gz
    # │   ├── b691deddf179ead0a006527330d86dde.sig.gz
    # │   ├── b7f087146f5cc3121477c29ff003e3d0.sig.gz
    # │   ├── c39c52d2d088348c950c2afe503b405b.sig.gz
    # │   ├── c9eb6a9d058df8036ad93bc45d5bf260.sig.gz
    # │   └── ce54d962851b0fdeefc624300036a133.sig.gz
    # ├── SOURMASH-MANIFEST.csv
    # ├── training_multisearch_result.csv
    # └── training_sig_files.txt
    # remove the files if they exist
    for f in expected_files:
        if exists(f):
            os.remove(f)
    #  python ../make_training_data_from_sketches.py --ref_file testdata/20_genomes_sketches.zip --ksize 31 --prefix 20_genomes_trained --outdir testdata/
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

    # then do the presence/absence estimation
    if exists(abundance_file):
        os.remove(abundance_file)
    # python ../run_YACHT.py --json testdata/20_genomes_trained_config.json --sample_file testdata/sample.sig.zip --out_file result.xlsx --outdir testdata/
    cmd = f"python {os.path.join(script_dir, 'run_YACHT.py')} --json {os.path.join(script_dir, 'gtdb_ani_thresh_0.95_config.json')} --sample_file {sample_sketches} --significance 0.99 --min_coverage 1 --outdir {data_dir} --out_filename {abundance_file}"
    # ~/pycharm/YACHT/tests/testdata$ tree 20_genomes_trained_intermediate_files/
    # 20_genomes_trained_intermediate_files/
    # ├── organism_sig_file.txt  # <-- new
    # ├── sample_multisearch_result.csv
    # ├── sample_sig_file.txt
    # ├── signatures
    # │   ├── 04212e93c2172d4df49dc5d8c2973d8b.sig.gz
    # │   ├── 04f2b0e94f2d1f1f5b8355114b70274e.sig.gz
    # │   ├── 0661ecab88c3d65d0f10e599a5ba1654.sig.gz
    # │   ├── 06ebe48d527882bfa9505aba8e31ae23.sig.gz
    # │   ├── 11fe9a00287c7ad086ebbc463724cf10.sig.gz
    # │   ├── 16c6c1d37259d83088ab3a4f5b691631.sig.gz
    # │   ├── 188d55801a78d4773cf6c0b46bca96ba.sig.gz
    # │   ├── 1a121dca600c6504e88252e81004f0cf.sig.gz
    # │   ├── 39ea7fd48ee7003587c9c763946d5d6e.sig.gz
    # │   ├── 45f2675c1dca4ef1a24a05f5b268adbb.sig.gz
    # │   ├── 7b312fffa3fb35440ba40203ba826c05.sig.gz
    # │   ├── 8fb9b1a69838a58cc4f31c1e42a5f189.sig.gz
    # │   ├── 92fb1b3e4baa6c474aff3efb84957687.sig.gz
    # │   ├── 96cb85214535b0f9723a6abc17097821.sig.gz
    # │   ├── a136145bee08846ed94c0406df3da2d4.sig.gz
    # │   ├── b691deddf179ead0a006527330d86dde.sig.gz
    # │   ├── b7f087146f5cc3121477c29ff003e3d0.sig.gz
    # │   ├── c39c52d2d088348c950c2afe503b405b.sig.gz
    # │   ├── c9eb6a9d058df8036ad93bc45d5bf260.sig.gz
    # │   └── ce54d962851b0fdeefc624300036a133.sig.gz
    # ├── SOURMASH-MANIFEST.csv
    # ├── training_multisearch_result.csv
    # └── training_sig_files.txt
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
