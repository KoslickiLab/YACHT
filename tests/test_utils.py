import subprocess
from os.path import exists
import os
import numpy as np
import pandas as pd
# add the parent directory to the path
import sys
cpath = os.path.dirname(os.path.realpath(__file__))
project_path = os.path.join(cpath,'..')
sys.path.append(project_path)
from yacht import utils
import sourmash
import unittest
import math
import json
import pytest
import tempfile
import gzip
import sys
import shutil

FILENAME = 'sample.sig.zip'
PATH_TO_YACHT_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/results')
PATH_TO_GENOME_TO_TAXID = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid.tsv')

def to_testing_data(file):
    return os.path.join(project_path, 'tests', os.path.join("testdata", file))


def test_load_signature_with_ksize1():
    # first, just try a *.sig file
    file = to_testing_data(FILENAME)
    sig = utils.load_signature_with_ksize(file, 31)
    # right type?
    assert type(sig) == sourmash.signature.FrozenSourmashSignature
    # can we do a simple operation on it?
    assert sig.jaccard(sig) == 1.0


def test_load_signature_with_ksize2():
    # wrong k-size
    file = to_testing_data(FILENAME)
    try:
        sig = utils.load_signature_with_ksize(file, 31)
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
    file = to_testing_data(FILENAME)
    sig = utils.load_signature_with_ksize(file, 31)
    sourmash.save_signatures([sig], open(to_testing_data('test.sig.zip'), 'wb'), compression=1)
    sig = utils.load_signature_with_ksize(to_testing_data('test.sig.zip'), 31)
    assert type(sig) == sourmash.signature.FrozenSourmashSignature
    assert sig.jaccard(sig) == 1.0

class TestGetColumnIndices(unittest.TestCase):
    def test_1(self):
        column_name_to_index = {
            "TAXID": 1,
            "RANK": 0,
            "PERCENTAGE": 2,
            "TAXPATH": 3,
            "TAXPATHSN": 4
        }
        indices = utils.get_column_indices(column_name_to_index)
        assert indices == (0, 1, 2, 3, 4)

    def test_2(self):
        column_name_to_index = {
            "RANK": 0,
            "PERCENTAGE": 2,
            "TAXPATH": 3,
            "TAXPATHSN": 4
        }
        with self.assertRaises(RuntimeError):
            utils.get_column_indices(column_name_to_index)

    def test_3(self):
        column_name_to_index = {
            "TAXID": 1,
            "PERCENTAGE": 2,
            "TAXPATH": 3,
            "TAXPATHSN": 4
        }
        with self.assertRaises(RuntimeError):
            utils.get_column_indices(column_name_to_index)

    def test_4(self):
        column_name_to_index = {
            "TAXID": 1,
            "RANK": 0,
            "TAXPATH": 3,
            "TAXPATHSN": 4
        }
        with self.assertRaises(RuntimeError):
            utils.get_column_indices(column_name_to_index)

    def test_5(self):
        column_name_to_index = {
            "TAXID": 1,
            "RANK": 0,
            "PERCENTAGE": 2,
            "TAXPATHSN": 4
        }
        with self.assertRaises(RuntimeError):
            utils.get_column_indices(column_name_to_index)

    def test_6(self):
        column_name_to_index = {
            "TAXID": 1,
            "RANK": 0,
            "PERCENTAGE": 2,
            "TAXPATH": 3
        }
        indices = utils.get_column_indices(column_name_to_index)
        assert indices[4] is None

class TestStandardizeOutput(unittest.TestCase):
    def test_everything_exists(self):

        yacht_output_dir = PATH_TO_YACHT_OUTPUT_DIR
        assert os.path.exists(yacht_output_dir)

        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)

        cmd = f"yacht convert --yacht_output_dir {yacht_output_dir} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert os.path.exists(os.path.join(outdir, 'cami_result.cami'))

    def test_wrong_yacht_output_dir(self):

        yacht_output_dir = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata_nonexisting')
        assert not os.path.exists(yacht_output_dir)

        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)

        cmd = f"yacht convert --yacht_output_dir {yacht_output_dir} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            _ = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_genome_to_taxid(self):

        yacht_output_dir = PATH_TO_YACHT_OUTPUT_DIR
        assert os.path.exists(yacht_output_dir)

        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid_nonexisting.tsv')
        assert not os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)

        cmd = f"yacht convert --yacht_output_dir {yacht_output_dir} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            _ = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_outdir(self):

        yacht_output_dir = PATH_TO_YACHT_OUTPUT_DIR
        assert os.path.exists(yacht_output_dir)

        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata_nonexisting')
        cmd = 'rm -rf ' + outdir
        try:
            subprocess.run(cmd, shell=True, check=True)
        except:
            pass
        assert not os.path.exists(outdir)

        cmd = f"yacht convert --yacht_output_dir {yacht_output_dir} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert os.path.exists(outdir)

        cmd = 'rm -rf ' + outdir
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0


if __name__ == '__main__':
    unittest.main()
