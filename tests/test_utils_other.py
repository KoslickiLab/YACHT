import unittest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
import sourmash
from yacht.utils import (load_signature_with_ksize, get_num_kmers,
                         get_info_from_single_sig,
                         get_column_indices,
                         create_output_folder)


class TestUtils(unittest.TestCase):

    def __init__(self, methodname: str = "runTest"):
        super().__init__(methodname)
        self.database = None

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_get_num_kmers(self):
        mean_abundance = 2.5
        hashes_len = 100
        scaled = 10000
        expected_num_kmers = int(mean_abundance * hashes_len * scaled)
        num_kmers = get_num_kmers(mean_abundance, hashes_len, scaled)
        self.assertEqual(num_kmers, expected_num_kmers)

    def test_get_info_from_single_sig(self):
        mock_sig = MagicMock()
        mock_sig.name = "test_signature"
        mock_sig.md5sum.return_value = "md5checksum"
        mock_sig.minhash.mean_abundance = 2.5
        mock_sig.minhash.hashes = {1: 1, 2: 1}
        mock_sig.minhash.scaled = 10000

        with patch('yacht.utils.load_signature_with_ksize', return_value=mock_sig):
            sig_info = get_info_from_single_sig("dummy_path", 31)[1:]
            self.assertEqual(sig_info, ("test_signature", "md5checksum", 2.5, 2, 10000))

    def test_get_column_indices(self):
        column_name_to_index = {
            "TAXID": 0,
            "RANK": 1,
            "PERCENTAGE": 2,
            "TAXPATH": 3,
            "TAXPATHSN": 4
        }
        expected_indices = (1, 0, 2, 3, 4)
        indices = get_column_indices(column_name_to_index)
        self.assertEqual(indices, expected_indices)

    def test_create_output_folder(self):
        outfolder = os.path.join(self.test_dir, 'outfolder')
        create_output_folder(outfolder)
        self.assertTrue(os.path.exists(outfolder))
