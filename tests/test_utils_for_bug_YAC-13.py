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
from  srcs.utils import load_signature_with_ksize, get_info_from_single_sig

def test_error_message_for_load_signature_with_ksize():
     # Reference to YAC-13 bug on sketching short sequence/genomes
     file = "tests/testdata/97_Silva_111_rep_set_euk_singleton.sig.zip"
     with pytest.raises(ValueError, match="Unable to calculate abundance mean. Please try sketch with '--scaled=1' or use alternative to sourmash."):
        load_signature_with_ksize(file, 31)
        
