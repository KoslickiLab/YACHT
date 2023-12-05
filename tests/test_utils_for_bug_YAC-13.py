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
import pytest
from  srcs.utils import remove_corr_organisms_from_ref, check_file_existence, get_cami_profile, get_column_indices, get_info_from_single_sig, collect_signature_info, run_multisearch


def test_error_message_for_get_info_from_single_sig():
     # Reference to YAC-13 bug on sketching short genomes
     file = "testdata/short_genomes_bug-YAC-13.sig.zip"
     with pytest.raises(ValueError, match="Empty sketch. Potential issues and suggestions: (1) The sketch is too big. Please try sketching with '--scaled=1', (2) Sequences are too small. Please use alternative to sourmash."):
         sig=get_info_from_single_sig(file, 31)

