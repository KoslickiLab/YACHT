import subprocess
from os.path import exists
import os
import numpy as np
import pandas as pd
# add the parent directory to the path
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from yacht import utils
import sourmash
import pytest
from  yacht.utils import load_signature_with_ksize

def test_error_message_for_load_signature_with_ksize():
     # Reference to YAC-13 bug on sketching short sequence/genomes
     file = "tests/testdata_bug_YAC13/extract_empty_hash.sig.zip"
     with pytest.raises(ValueError, match=f"Empty sketch in signature. This may be due to too high of a scale factor, please reduce it, eg. --scaled=1, and try again."):
        load_signature_with_ksize(file, 31)

        
