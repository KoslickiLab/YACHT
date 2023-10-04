import subprocess
from os.path import exists
import os
import numpy as np
import pandas as pd
# add the parent directory to the path
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from srcs import utils
import sourmash


def to_testing_data(file):
    return os.path.join('tests', os.path.join("testdata", file))


def test_load_hashes_to_index():
     # the *hash_to_col_idx.pkl files contain a pickle list of key value pairs with keys the hash values and values
     # the index of the row that they appear in the npz matrix
     file = to_testing_data("integration_test_hash_to_col_idx.pkl")
     hashes = utils.load_hashes_to_index(file)
     assert type(hashes) == dict
     assert len(hashes) == 63888
     assert np.allclose(np.sort(list(hashes.values())), range(0, len(hashes)))


def test_load_signature_with_ksize1():
    # first, just try a *.sig file
    file = to_testing_data("sample.sig")
    sig = utils.load_signature_with_ksize(file, 31)
    # right type?
    assert type(sig) == sourmash.signature.FrozenSourmashSignature
    # can we do a simple operation on it?
    assert sig.jaccard(sig) == 1.0


def test_load_signature_with_ksize2():
    # wrong k-size
    file = to_testing_data("sample.sig")
    try:
        sig = utils.load_signature_with_ksize(file, 21)
    except ValueError:
        pass
    # wrong file type
    file = to_testing_data("foobar")
    try:
        sig = utils.load_signature_with_ksize(file, 31)
    except ValueError:
        pass
    # too many files
    try:
        sig = utils.load_signature_with_ksize(to_testing_data("20_genomes_sketches.zip"), 31)
    except ValueError:
        pass


def test_load_signature_with_ksize3():
    # different kind of format
    file = to_testing_data("sample.sig")
    sig = utils.load_signature_with_ksize(file, 31)
    sourmash.save_signatures([sig], open(to_testing_data('test.sig.zip'), 'wb'), compression=1)
    sig = utils.load_signature_with_ksize(to_testing_data('test.sig.zip'), 31)
    assert type(sig) == sourmash.signature.FrozenSourmashSignature
    assert sig.jaccard(sig) == 1.0







