import unittest
from unittest.mock import patch

import requests

import yacht.download_pretrained_ref_db as download_script
import argparse
import os


class TestDownloadPretrainedRefDb(unittest.TestCase):
    GENBANK_VERSION = 'genbank-2022.03'
    PRETRAINED_ZIP = 'genbank-2022.03-bacteria-k31_0.95_pretrained.zip'

    def setUp(self):
        self.parser = argparse.ArgumentParser()
        download_script.add_arguments(self.parser)

    def test_add_arguments(self):
        args = self.parser.parse_args(['--database', 'genbank', '--db_version', self.GENBANK_VERSION,
                                       '--ani_thresh', '0.95', '--outfolder', '.'])
        self.assertEqual(args.database, 'genbank')
        self.assertEqual(args.db_version, self.GENBANK_VERSION)
        self.assertEqual(args.ani_thresh, 0.95)
        self.assertEqual(args.outfolder, '.')

    @patch('yacht.download_pretrained_ref_db.requests.get')
    def test_fetch_zenodo_records_success(self, mock_get):
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {'hits': {'hits': [{'title': 'test-db'}]}}
        result = download_script.fetch_zenodo_records()
        self.assertEqual(result, [{'title': 'test-db'}])

    @patch('yacht.download_pretrained_ref_db.requests.get')
    def test_fetch_zenodo_records_failure(self, mock_get):
        mock_get.side_effect = requests.exceptions.RequestException('Failed to fetch')
        result = download_script.fetch_zenodo_records()
        self.assertEqual(result, [])

    @patch('yacht.download_pretrained_ref_db.requests.get')
    @patch('yacht.download_pretrained_ref_db.download_file')
    @patch('yacht.download_pretrained_ref_db.fetch_zenodo_records')
    @patch('yacht.download_pretrained_ref_db.unzip_file')
    def test_main(self, mock_unzip, mock_fetch_records, mock_download_file, mock_get):
        args = argparse.Namespace(database='genbank', db_version=self.GENBANK_VERSION, ncbi_organism="bacteria",
                                  ani_thresh=0.95, k=31, outfolder='.')
        mock_fetch_records.return_value = [{'title': self.PRETRAINED_ZIP, 'files': [
            {'key': self.PRETRAINED_ZIP, 'links': {'self': 'https://zenodo.url'}}]}]
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {'hits': {'hits': [{'title': 'test-db'}]}}
        mock_download_file.return_value = True
        download_script.main(args)
        mock_unzip.assert_called_with(os.path.join('.', self.PRETRAINED_ZIP), '.')
