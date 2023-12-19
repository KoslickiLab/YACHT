import unittest
import argparse
import os
import tempfile
import shutil
import zipfile
from unittest.mock import patch, MagicMock
import pandas as pd

from yacht import make_training_data_from_sketches


class TestMakeTrainingDataFromSketches(unittest.TestCase):
    TEST_ZIP_FILE_NAME = 'test_ref_file.zip'

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.test_zip_file = os.path.join(self.temp_dir, self.TEST_ZIP_FILE_NAME)

        with zipfile.ZipFile(self.test_zip_file, 'w') as dummy_zip:
            dummy_zip.writestr('SOURMASH-MANIFEST.csv', 'test,content\n')
            dummy_zip.writestr('signatures/test.sig.gz', 'dummy signature content')

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_main_valid_input(self):
        test_args = argparse.Namespace(
            ref_file=self.test_zip_file,
            ksize=21,
            num_threads=4,
            ani_thresh=0.95,
            prefix='test_prefix',
            outdir=self.temp_dir,
            force=False
        )

        with patch('yacht.make_training_data_from_sketches.utils') as mock_utils:
            mock_utils.check_file_existence = MagicMock()
            mock_utils.collect_signature_info = MagicMock(return_value={'sig1': ['sig1', 1]})
            mock_utils.run_multisearch = MagicMock(return_value=pd.DataFrame())
            mock_df = MagicMock(spec=pd.DataFrame)
            mock_utils.remove_corr_organisms_from_ref = MagicMock(return_value=(mock_df, mock_df))

            make_training_data_from_sketches.main(test_args)

            ref_file_error_msg = f'Reference database zip file {self.test_zip_file} does not exist.'
            mock_utils.check_file_existence.assert_called_with(self.test_zip_file, ref_file_error_msg)

    def test_main_invalid_zip_file(self):
        test_args = argparse.Namespace(
            ref_file=os.path.join(self.temp_dir, 'invalid_file.txt'),
            ksize=21,
            num_threads=4,
            ani_thresh=0.95,
            prefix='test_prefix',
            outdir=self.temp_dir,
            force=False
        )

        with self.assertRaises(ValueError):
            make_training_data_from_sketches.main(test_args)


if __name__ == '__main__':
    unittest.main()
