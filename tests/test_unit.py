# Unit tests
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from hypothesis_recovery import get_alt_mut_rate


def test_get_alt_mut_rate_1():
    assert get_alt_mut_rate(100, 10000, 21, significance=0.99) == -1