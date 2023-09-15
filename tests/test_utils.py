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
    return os.path.join("testdata", file)


def test_load_hashes():
    # the *hash_to_col_idx.pkl files contain a pickle list of key value pairs with keys the hash values and values
    # the index of the row that they appear in the npz matrix
    file = to_testing_data("integration_test_hash_to_col_idx.pkl")
    hashes = utils.load_hashes_to_index(file)
    assert type(hashes) == dict
    assert len(hashes) == 63888
    assert np.allclose(np.sort(list(hashes.values())), range(0, len(hashes)))


def test_load_signature_with_ksize():
    # first, just try a *.sig file
    file = to_testing_data("sample.sig")
    sig = utils.load_signature_with_ksize(file, 31)
    assert type(sig) == sourmash.signature.FrozenSourmashSignature




