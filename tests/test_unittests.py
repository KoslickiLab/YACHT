import json
import math
import pytest
import pandas as pd
import os
import tempfile
import gzip
import sys
import shutil
cpath = os.path.dirname(os.path.realpath(__file__))
project_path = os.path.join(cpath,'..')
sys.path.append(project_path)
from yacht.hypothesis_recovery_src import single_hyp_test,  get_alt_mut_rate
from yacht.utils import check_file_existence, get_cami_profile, get_column_indices, get_info_from_single_sig, collect_signature_info

@pytest.fixture
def test_output_files():
    filename = 'test_output'
    yield filename
    if os.path.exists(filename):
        os.remove(filename)

tmp_dir = f"{project_path}/tests/unittests_data/test_tmp"

hashes_data = {'hash1': 1, 'hash2': 2, 'hash3': 3}
ksize = 31

def test_check_file_existence():
    dont_exist = 'File does not exist'
    existing_file = os.path.join(tmp_dir, "existing_file.txt")
    with open(existing_file, 'w') as f:
        f.write("Test content")
    assert check_file_existence(existing_file, dont_exist) is None

    non_existing_file = os.path.join(tmp_dir, "non_existing_file.txt")
    with pytest.raises(ValueError, match=dont_exist):
        check_file_existence(non_existing_file, dont_exist)
        
def test_get_column_indices():
    column_name_to_index = {
        "TAXID": 1,
        "RANK": 0,
        "PERCENTAGE": 2,
        "TAXPATH": 3,
        "TAXPATHSN": 4
    }
    indices = get_column_indices(column_name_to_index)
    assert indices == (0, 1, 2, 3, 4)
    
def test_get_cami_profile():
    file_path = os.path.join(project_path, 'tests/testdata/sample_cami.txt')
    with open(file_path, 'r') as file:
        sample_cami_content = file.readlines()
    
    profiles = get_cami_profile(sample_cami_content)

    expected_header = {
        'SAMPLEID': 'CAMI_LOW_S001', 
        'VERSION': '0.9.1', 
        'RANKS': 'superkingdom|phylum|class|order|family|genus|species|strain', 
        'TAXONOMYID': 'ncbi-taxonomy_DATE', 
        '__PROGRAM__': 'unknown'
    }

    assert len(profiles) == 1
    sample_id, header, profile = profiles[0]

    assert sample_id == "CAMI_LOW_S001"
    assert header == expected_header
    assert len(profile) == 23

    prediction1 = profile[0]
    assert prediction1.rank == 'superkingdom'
    assert prediction1.taxid == '2'
    assert math.isclose(prediction1.percentage, 29.183763, abs_tol=1e-6)
    assert prediction1.taxpath == '2'
    assert prediction1.taxpathsn == 'Bacteria'

    prediction2 = profile[1]
    assert prediction2.rank == 'phylum'
    assert prediction2.taxid == '201174'
    assert math.isclose(prediction2.percentage, 4.638241, rel_tol=1e-6)
    assert prediction2.taxpath == '2|201174'
    assert prediction2.taxpathsn == 'Bacteria|Actinobacteria'

def test_get_alt_mut_rate():
    nu = 10
    thresh = 5
    ksize = 31
    significance = 0.99
    result = get_alt_mut_rate(nu, thresh, ksize, significance)
    expected_result = 0.047902071844405425
    assert math.isclose(result, expected_result, rel_tol=1e-6, abs_tol=1e-6)

def test_get_alt_mut_rate_zero_nu():
    nu = 0
    thresh = 5
    ksize = 31
    significance = 0.99
    result = get_alt_mut_rate(nu, thresh, ksize, significance)
    expected_result = -1
    assert result == expected_result

def test_get_alt_mut_rate_large_thresh():
    nu = 10
    thresh = 20
    ksize = 31
    significance = 0.99
    result = get_alt_mut_rate(nu, thresh, ksize, significance)
    expected_result = -1
    assert result == expected_result

def test_get_info_from_single_sig():
    sig_list_file = f'{project_path}/gtdb_ani_thresh_0.95_intermediate_files/training_sig_files.tsv'
    with open(sig_list_file, 'r') as file:
        lines = file.readlines()
        if lines:
            sig_file_path = [line.strip() for line in lines if "96cb85214535b0f9723a6abc17097821.sig.gz" in line][0]
        else:
            raise IOError("Signature list file is empty")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_sig_file = os.path.join(tmpdir, os.path.basename(sig_file_path))

        with gzip.open(sig_file_path, 'rb') as f_in:
            with open(tmp_sig_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        ksize = 0
        result = get_info_from_single_sig(tmp_sig_file, ksize)[1:]

        expected_name = "VIKJ01000003.1 Chitinophagaceae bacterium isolate X1_MetaBAT.39 scaffold_1008, whole genome shotgun sequence"
        expected_md5sum = "96cb85214535b0f9723a6abc17097821"
        expected_mean_abundance = 1.0
        expected_hashes_len = 1984
        expected_scaled = 1000

        assert result[0] == expected_name
        assert result[1] == expected_md5sum
        assert abs(result[2] - expected_mean_abundance) < 0.01
        assert result[3] == expected_hashes_len
        assert result[4] == expected_scaled

def test_collect_signature_info():
    num_threads = 2
    ksize = 0
    path_to_temp_dir = f'{project_path}/gtdb_ani_thresh_0.95_intermediate_files/' 
    result = collect_signature_info(num_threads, ksize, path_to_temp_dir)

    with open(f'{project_path}/tests/unittests_data/test_collect_signature_info_data.json', 'r') as file:
        expectations = json.load(file)

    for expectation in expectations.keys():
        assert expectation in result
        actual_info = result[expectation]
        assert expectations[expectation] == list(actual_info[:-1])

def test_single_hyp_test():
    exclusive_hashes_info_org = (100, 90)
    ksize = 31

    result = single_hyp_test(exclusive_hashes_info_org, ksize)

    in_sample_est, p_val, num_exclusive_kmers, num_exclusive_kmers_coverage, num_matches, \
    acceptance_threshold_with_coverage, actual_confidence_with_coverage, alt_confidence_mut_rate_with_coverage = result

    assert isinstance(in_sample_est, int)
    assert isinstance(p_val, float)
    assert isinstance(num_exclusive_kmers, int)
    assert isinstance(num_exclusive_kmers_coverage, int)
    assert isinstance(num_matches, int)
    assert isinstance(acceptance_threshold_with_coverage, float)
    assert isinstance(actual_confidence_with_coverage, float)
    assert isinstance(alt_confidence_mut_rate_with_coverage, float)


if __name__ == '__main__':
    pytest.main()
