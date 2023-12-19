# Unit tests
import numpy as np
import sys
import os
cpath = os.path.dirname(os.path.realpath(__file__))
project_path = os.path.join(cpath,'..')
sys.path.append(project_path)
from yacht.hypothesis_recovery_src import get_alt_mut_rate


def test_get_alt_mut_rate_1():
    # All these are calculated with Steve's code
    assert get_alt_mut_rate(100, 10000, 21, significance=0.99) == -1
    assert np.isclose(get_alt_mut_rate(10, 0, 21), 0.28015945851802826)
    assert np.isclose(get_alt_mut_rate(10, 0, 31), 0.19963312102481723)
    assert np.isclose(get_alt_mut_rate(10, 5, 21), 0.0698992155957967)
    assert np.isclose(get_alt_mut_rate(10, 5, 31), 0.047902071848511696)
    assert np.isclose(get_alt_mut_rate(10, 9, 21), 0.02169068099465221)
    assert np.isclose(get_alt_mut_rate(100, 10, 11), 0.2397729973308742)
    assert np.isclose(get_alt_mut_rate(1000, 0, 1), 0.9999899497147453)
