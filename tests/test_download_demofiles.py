import unittest
from unittest.mock import patch, mock_open
import yacht.download_demofiles as demo_files
import argparse
import os

class TestDownloadDemoFiles(unittest.TestCase):
    REAL_DEMO_URL = "https://raw.githubusercontent.com/KoslickiLab/YACHT/main/demo"
    FILE_LIST_URL = "https://api.github.com/repos/KoslickiLab/YACHT/contents/demo"
    FILE_NAME = 'file1.txt'

    def setUp(self):
        self.parser = argparse.ArgumentParser()
        demo_files.add_arguments(self.parser)

    def test_add_arguments(self):
        args = self.parser.parse_args(['--outfolder', 'demo'])
        self.assertEqual(args.outfolder, 'demo')

    @patch('yacht.download_demofiles.requests.get')
    def test_fetch_file_list_from_github_success(self, mock_get):
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = [{'path': os.path.join('demo', self.FILE_NAME), 'type': 'file'}]
        result = demo_files.fetch_file_list_from_github()
        self.assertEqual(result, [self.FILE_NAME])

    @patch('yacht.download_demofiles.fetch_file_list_from_github')
    @patch('yacht.download_demofiles.download_file')
    def test_download_demo_files(self, mock_download_file, mock_fetch_list):
        mock_fetch_list.return_value = [self.FILE_NAME]
        demo_files.download_demo_files('demo')
        expected_url = os.path.join(self.REAL_DEMO_URL, self.FILE_NAME)
        mock_download_file.assert_called_with(expected_url, os.path.join('demo', self.FILE_NAME))

    @patch('yacht.download_demofiles.requests.get')
    def test_fetch_file_list_from_github_failure(self, mock_get):
        mock_get.side_effect = Exception('Failed to fetch')
        with self.assertRaises(Exception):
            demo_files.fetch_file_list_from_github()

    @patch('yacht.download_demofiles.requests.get')
    @patch('builtins.open', new_callable=mock_open)
    def test_download_file(self, mock_file_open, mock_get):
        mock_get.return_value.status_code = 200
        mock_get.return_value.content = b'data'
        demo_files.download_file(self.REAL_DEMO_URL, 'path/to/file')
        mock_file_open.assert_called_with('path/to/file', 'wb')

    @patch('yacht.download_demofiles.download_demo_files')
    def test_main(self, mock_download_demo_files):
        args = argparse.Namespace(outfolder='demo')
        demo_files.main(args)
        mock_download_demo_files.assert_called_with('demo')
