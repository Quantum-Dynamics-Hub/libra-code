"""
Unit and regression test for the workflows/nbra module in the Libra package
"""

from libra_py.workflows.nbra import step3_auxiliary
import pytest
import sys
import os
 
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



@pytest.mark.parametrize("data_dim, gs_index, excitation_type, expected_result", [
    ( 4, 0, 0, [[1,-1], [1,-2]] ),   
    ( 4, 0, 1, [[1,-1], [-1,2]] ),
    ( 4, 0, 2, [[1,-1], [1,-2], [-1,2]] ),
    ( 4, 1, 0, [[1,-1,2,-2]]),
    ( 8, 0, 0, [[1,-1], [1,-2], [1,-3], [1,-4]] ),
    ( 8, 1, 0, [[1,-1,2,-2], [1,2,-2,-3], [1,2,-2,-4], [1,-1,2,-3], [1,-1,2,-4]] ),
    ( 8, 2, 0, [[1,-1,2,-2,3,-3], [1,2,-2,3,-3,-4], [1,-1,2,3,-3,-4], [1,-1,2,-2,3,-4]] ),
    ( 8, 3, 0, [[1,-1,2,-2,3,-3,4,-4]] ),
    ( 8, 0, 2, [[1,-1], [1,-2], [1,-3], [1,-4], [-1,2], [-1,3], [-1,4]] ),
    ( 8, 3, 2, [[1,-1,2,-2,3,-3,4,-4]] ),
])

def test_build_SD_pyxaid(data_dim, gs_index, excitation_type, expected_result):
    """Tests that the build_SD_pyxaid function builds the SD basis we expect."""   
    computed_result = step3_auxiliary.build_SD_pyxaid(data_dim, gs_index, excitation_type)
    assert expected_result == computed_result




@pytest.mark.parametrize("data_dim, cbm_alpha_index, cbm_beta_index, excitation_type, expected_result", [
    ( 6, 2, 5, 1, [ [1,2,-4,-5], [1,2,-6,-5], [1,2,-4,-6] ] ),
    ( 6, 2, 4, 1, [ [1,2,-4], [1,2,-5], [1,2,-6] ] ),
    ( 6, 1, 4, 0, [ [1,-4], [2,-4], [3,-4] ] ),
    ( 4, 1, 3, 2, [ [1,-3], [2,-3], [1,-4] ] ),
    ( 6, 2, 4, 2, [ [1,2,-4], [3,2,-4], [1,3,-4], [1,2,-5], [1,2,-6] ] ),
    (10, 2, 7, 0, [ [1,2,-6,-7], [3,2,-6,-7], [4,2,-6,-7], [5,2,-6,-7], [1,3,-6,-7], [1,4,-6,-7], [1,5,-6,-7] ] ),
    ( 4, 2, 4, 2, [ [1,2,-3,-4] ] ),
    ( 4, 2, 4, 1, [ [1,2,-3,-4] ] ),
    ( 4, 2, 4, 0, [ [1,2,-3,-4] ] ),
])
def test_build_SD(data_dim, cbm_alpha_index, cbm_beta_index, excitation_type, expected_result):
    """Tests that the build_SD_pyxaid function builds the SD basis we expect."""
    computed_result = step3_auxiliary.build_SD(data_dim, cbm_alpha_index, cbm_beta_index, excitation_type)
    assert expected_result == computed_result



