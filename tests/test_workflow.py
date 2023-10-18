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
    # In testdata/
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
        assert os.stat(f).st_size > 400
    # then do the presence/absence estimation
    if exists(abundance_file):
        os.remove(abundance_file)
    # python ../run_YACHT.py --json testdata/20_genomes_trained_config.json --sample_file testdata/sample.sig.zip --out_file result.xlsx
    cmd = f"python {os.path.join(script_dir, 'run_YACHT.py')} --json {os.path.join(data_dir, '20_genomes_trained_config.json')} --sample_file {sample_sketches} --significance 0.99 --min_coverage 0.001 --out_file {os.path.join(data_dir,abundance_file)} --show_all"
    print(cmd)
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
    res = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; sourmash sketch fromfile ref_paths.csv -p dna,k=31,scaled=1000,abund -o ref.sig.zip --force-output-already-exists"
    res = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../make_training_data_from_sketches.py --force --ref_file ref.sig.zip --ksize 31 --num_threads 1 --ani_thresh 0.95 --prefix 'demo_ani_thresh_0.95' --outdir ./"
    res = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../run_YACHT.py --json demo_ani_thresh_0.95_config.json --sample_file sample.sig.zip --significance 0.99 --num_threads 1 --min_coverage_list 1 0.6 0.2 0.1 --out_filename result.xlsx"
    res = subprocess.run(cmd, shell=True, check=True)
    cmd = "cd demo; python ../srcs/standardize_yacht_output.py --yacht_output result.xlsx --sheet_name min_coverage0.2 --genome_to_taxid toy_genome_to_taxid.tsv --mode cami --sample_name 'MySample' --outfile_prefix cami_result --outdir ./"
    res = subprocess.run(cmd, shell=True, check=True)


