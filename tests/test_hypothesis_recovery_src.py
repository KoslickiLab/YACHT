import unittest
from unittest.mock import patch, MagicMock, mock_open
import pandas as pd
import sourmash

from yacht.hypothesis_recovery_src import (
    get_organisms_with_nonzero_overlap,
    get_exclusive_hashes,
    single_hyp_test,
    hypothesis_recovery
)


class TestHypothesisRecoverySrc(unittest.TestCase):

    def setUp(self):
        self.mock_manifest = pd.DataFrame({
            'organism_name': ['org1', 'org2'],
            'md5sum': ['md5_1', 'md5_2'],
            'num_unique_kmers_in_genome_sketch': [1000, 2000],
            'num_total_kmers_in_genome_sketch': [5000, 6000],
            'genome_scale_factor': [0.1, 0.2],
            'num_exclusive_kmers_in_sample_sketch': [100, 200],
            'num_total_kmers_in_sample_sketch': [700, 800],
            'sample_scale_factor': [0.3, 0.4],
            'min_coverage': [0.5, 0.6]
        })
        self.sample_sig = MagicMock(spec=sourmash.SourmashSignature)
        self.sample_info_set = ('sample_file_path', self.sample_sig)

    @patch('yacht.hypothesis_recovery_src.get_organisms_with_nonzero_overlap')
    @patch('yacht.hypothesis_recovery_src.get_exclusive_hashes')
    def test_hypothesis_recovery(self, mock_get_exclusive_hashes, mock_get_organisms_with_nonzero_overlap):
        mock_get_organisms_with_nonzero_overlap.return_value = ['org1', 'org2']
        mock_get_exclusive_hashes.return_value = ([(10, 5)], self.mock_manifest)

        result = hypothesis_recovery(self.mock_manifest, self.sample_info_set,
                                     'path_to_genome_temp_dir', [0.1, 0.2], 1, 31)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    @patch('pandas.read_csv')
    @patch('yacht.hypothesis_recovery_src.os.system')
    @patch('yacht.hypothesis_recovery_src.os.path.join')
    @patch('yacht.hypothesis_recovery_src.os.listdir')
    @patch('builtins.open', new_callable=mock_open, create=True)
    @patch('yacht.hypothesis_recovery_src.zipfile.ZipFile')
    def test_get_organisms_with_nonzero_overlap(self, mock_zipfile, _, mock_os_listdir,
                                                mock_os_path_join, mock_os_system, mock_read_csv):
        mock_os_listdir.return_value = ['sig_file']
        mock_os_path_join.return_value = 'joined_path'
        mock_os_system.return_value = 0
        mock_read_csv.return_value = pd.DataFrame({'match_name': ['org1', 'org2']})

        mock_zip_file_instance = mock_zipfile.return_value.__enter__.return_value
        mock_zip_file_instance.extractall.return_value = None

        result = get_organisms_with_nonzero_overlap(self.mock_manifest, 'sample_file.zip', 10, 31, 4,
                                                    '/path/to/genome_temp_dir', '/path/to/sample_temp_dir')

        self.assertIsInstance(result, list)

    @patch('yacht.hypothesis_recovery_src.load_signature_with_ksize')
    def test_get_exclusive_hashes(self, mock_load_signature_with_ksize):
        mock_load_signature_with_ksize.return_value = MagicMock(spec=sourmash.SourmashSignature)

        result = get_exclusive_hashes(self.mock_manifest, ['org1', 'org2'], self.sample_sig, 31,
                                      '/path/to/genome_temp_dir')

        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result[0], list)
        self.assertIsInstance(result[1], pd.DataFrame)

    def test_single_hyp_test(self):
        result = single_hyp_test((10, 5), 31, 0.99, 0.95)

        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 8)
