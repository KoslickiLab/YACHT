import unittest
from unittest.mock import patch, mock_open

import requests

import yacht.download_default_ref_db
import argparse
import os


class TestDownloadScript(unittest.TestCase):
    URL = "https://some.url"
    FILE_PATH = "path/to/file"
    FOLDER_PATH = "path/to/folder"

    def setUp(self):
        self.parser = argparse.ArgumentParser()
        yacht.download_default_ref_db.add_arguments(self.parser)

    def test_add_arguments(self):
        args = self.parser.parse_args(['--database', 'genbank', '--db_version', 'genbank-2022.03'])
        self.assertEqual(args.database, 'genbank')
        self.assertEqual(args.db_version, 'genbank-2022.03')

    def test_generate_download_url_genbank(self):
        args = argparse.Namespace(database="genbank", db_version="genbank-2022.03", ncbi_organism="virus", k=31)
        url = yacht.download_default_ref_db.generate_download_url(args)
        expected_url = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2022.03/genbank-2022.03-viral-k31.zip"
        self.assertEqual(url, expected_url)

    def test_invalid_genbank_version(self):
        invalid_version = "genbank-2023.99"  # Example invalid version
        args = argparse.Namespace(database="genbank", db_version=invalid_version, ncbi_organism="virus", k=31)
        url = yacht.download_default_ref_db.generate_download_url(args)
        self.assertIsNone(url)

    @patch('yacht.download_default_ref_db.requests.get')
    def test_download_file_success(self, mock_get):
        mock_get.return_value.status_code = 200
        mock_get.return_value.content = b"file content"
        with patch('builtins.open', mock_open()) as mock_file:
            result = yacht.download_default_ref_db.download_file(self.URL, self.FILE_PATH)
            mock_file.assert_called_with(self.FILE_PATH, 'wb')
            self.assertTrue(result)

    @patch('yacht.download_default_ref_db.requests.get')
    def test_download_file_failure(self, mock_get):
        mock_get.side_effect = requests.exceptions.RequestException('Failed to download')
        result = yacht.download_default_ref_db.download_file(self.URL, self.FILE_PATH)
        self.assertFalse(result)

    @patch('yacht.download_default_ref_db.download_file')
    @patch('yacht.download_default_ref_db.generate_download_url')
    @patch('yacht.download_default_ref_db.check_download_args')
    @patch('yacht.download_default_ref_db.create_output_folder')
    def test_main(self, mock_create_folder, mock_check_args, mock_generate_url, mock_download_file):
        mock_generate_url.return_value = self.URL
        mock_download_file.return_value = True
        args = argparse.Namespace(outfolder=self.FOLDER_PATH, database="genbank", db_version="genbank-2022.03")
        yacht.download_default_ref_db.main(args)
        mock_check_args.assert_called_once()
        mock_create_folder.assert_called_with(self.FOLDER_PATH)
        mock_generate_url.assert_called_with(args)
        mock_download_file.assert_called_with(self.URL, os.path.join(self.FOLDER_PATH, os.path.basename(self.URL)))
