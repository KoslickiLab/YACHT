import subprocess, os
import tempfile
import json
from os.path import exists

def make_train_fasta():
    fasta_content = [
        ">Genome1",
        "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        ">Genome2",
        "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
        ">Genome3",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
    ]

    fasta_filename = "example.fasta"

    with open(fasta_filename, "w") as fasta_file:
        for line in fasta_content:
            fasta_file.write(line + "\n")

def test_sourmash_sketch_command():
    with tempfile.TemporaryDirectory() as tmp_dir:
        make_train_fasta()

        fasta_file = "example.fasta"
        output_file = os.path.join(tmp_dir, "training_database.sig.zip")
        cmd = [
            "sourmash", "sketch", "dna", "-f", "-p", "k=31,scaled=1000,abund", "--singleton", fasta_file, "-o", output_file
        ]

        subprocess.run(cmd, check=True)

        assert os.path.isfile(output_file)

def test_make_training_data_from_sketches():
    ref_file = 'tests/testdata/20_genomes_sketches.zip'
    ksize = '31'
    ani_thresh = '0.95'
    prefix = 'gtdb_ani_thresh_0.95'
    config_file = f'{prefix}_config.json'
    processed_manifest_file = f'{prefix}_processed_manifest.tsv'
    intermediate_files_dir = f'{prefix}_intermediate_files'

    command = [
        'yacht', 'train',
        '--ref_file', ref_file,
        '--ksize', ksize,
        '--prefix', prefix,
        '--ani_thresh', ani_thresh,
        '--outdir', './',
        '--force',
    ]

    subprocess.run(command)

    assert os.path.isfile(config_file)
    assert os.path.isfile(processed_manifest_file)
    assert os.path.isdir(intermediate_files_dir)

    with open(config_file, 'r') as f:
        config = json.load(f)
        assert config['ksize'] == int(ksize)
        assert config['ani_thresh'] == float(ani_thresh)

def test_run_yacht():
    cmd = "yacht run --json gtdb_ani_thresh_0.95_config.json --sample_file 'tests/testdata/sample.sig.zip' --significance 0.99 --min_coverage_list 1 0.6 0.2 0.1"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0

    assert exists('result.xlsx')

def test_run_train_sig():
    ksize = '31'
    ani_thresh = '0.95'
    prefix = 'gtdb_ani_thresh_0.95'
    config_file = f'{prefix}_config.json'
    processed_manifest_file = f'{prefix}_processed_manifest.tsv'
    intermediate_files_dir = f'{prefix}_intermediate_files'
    cmd = "yacht train --ref_file 'tests/testdata/sample.sig' --ksize 31 --prefix 'gtdb_ani_thresh_0.95' --ani_thresh 0.95 --outdir ./ --force"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0

    assert os.path.isfile(config_file)
    assert os.path.isfile(processed_manifest_file)
    assert os.path.isdir(intermediate_files_dir)

    with open(config_file, 'r') as f:
        config = json.load(f)
        assert config['ksize'] == int(ksize)
        assert config['ani_thresh'] == float(ani_thresh)
        
def test_run_pretrained_ref_db():
    cmd = "yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 31 --ani_thresh 0.9995 --outfolder ./"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0

def test_run_yacht_pretrained_ref_db():
    cmd = "yacht run --json ./gtdb-rs214-reps.k31_0.9995_pretrained/gtdb-rs214-reps.k31_0.9995_config.json --sample_file tests/testdata/sample.sig.zip --significance 0.99 --num_threads 32 --min_coverage_list 1 0.6 0.2 0.1 --out ./result_pretrained.xlsx"
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0

    assert exists('result_pretrained.xlsx')

