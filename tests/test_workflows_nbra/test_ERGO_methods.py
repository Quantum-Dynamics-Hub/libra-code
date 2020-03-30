"""
Unit and regression test for the workflows/nbra module and its dependencies in the Libra package.
"""

from libra_py import ERGO_methods
import pytest
import sys
import os
 
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



@pytest.mark.parametrize("expected_result", [
    ( [[-0.275036, -0.272319, -0.250507, -0.240239, -0.22841, -0.0277049, -0.00211619, -0.00014508, 0.000937045, 0.0157079]], [5, 5]  ),   
])

def test_read_spectrum_restricted(expected_result):
    """Tests that the read_spectrum_restricted function reads ErgoSCF output files as we expect
       The pytest folder should have the following files:
           1. occupied_spectrum.txt
           2. unoccupied_spectrum.txt  
    """   
    computed_result = ERGO_methods.read_spectrum_restricted()
    assert expected_result == computed_result


