"""
Unit and regression test for the workflows/nbra module in the Libra package
"""

from libra_py.workflows.nbra import mapping
import pytest
import sys
import os
 
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *




def get_data(filename, data_dim):
    X = MATRIX(data_dim, data_dim)
    X.Load_Matrix_From_File(filename)
    return X




@pytest.mark.parametrize("inp, nbasis, do_sort, expected_result", [
    ( [1,2,-5,-6], 8, True, [0,1,4,5] ),   
    ( [1,3,-5,-6], 8, True, [0,2,4,5] ),
    ( [3,1,-5,-6], 8, True, [0,2,4,5] ),
    ( [1,2,3,-5],  8, True, [0,1,2,4] ),
    ( [-5,-6,-7,1],8, True, [0,4,5,6] ),
])
def test_sd2indx(inp, nbasis, do_sort, expected_result):
    """Tests that the sd2indx function indexes our input SDs in accordance to the step2 data format"""   
    computed_result = mapping.sd2indx(inp, nbasis, do_sort)
    assert expected_result == computed_result




@pytest.mark.parametrize("inp, nbasis, do_sort", [
    ( [1,2,10,-5], 8, True ),
    ( [8,9,-5,-6], 8, True ),
])
def test_sd2indx_exit(inp, nbasis, do_sort):
    """Tests that the sd2indx function sucesfully extis upon a rightful error"""
    with pytest.raises(SystemExit):
        mapping.sd2indx(inp, nbasis, do_sort)




@pytest.mark.parametrize("SD1, SD2, S, expected_result", [  
    ( [1,2,-5,-6], [1,2,-5,-6], get_data("supporting_files/overlap_matrix_example.txt", 8), 1.0+0.0j ),
    ( [1,2,-5,-6], [3,4,-7,-8], get_data("supporting_files/overlap_matrix_example.txt", 8), 0.0+0.0j ),
])
def test_ovlp_arb(SD1, SD2, S, expected_result):
    """Tests that the ovlp_arb function correctly computes the overlap of SD in accordance to the example kohn-sham data in the file overlap_matrix_example"""
    computed_result = mapping.ovlp_arb(SD1, SD2, S)
    assert expected_result == round(computed_result.real, 6) + round(computed_result.imag, 6) * 1j 




@pytest.mark.parametrize("SD1, Hvib_re, expected_result", [
    ( [1,2,-5,-6], get_data("supporting_files/hvib_matrix_real_example.txt", 8), 1.38909+0.0j ),
    ( [1,2,-5,-8], get_data("supporting_files/hvib_matrix_real_example.txt", 8), 1.43555+0.0j ),
    ( [1,2,3,4],   get_data("supporting_files/hvib_matrix_real_example.txt", 8), 1.48560+0.0j ),
])
def test_energy_arb(SD1, Hvib_re, expected_result):
    """Tests that the energy_arb function correctly computes the energy of SD in accordance to the example kohn-sham data in the file hvib_matrix_real_example
       This energy should be the sum of the energies of the individual kohn-sham spin-orbitals
    """
    computed_result = mapping.energy_arb(SD1, Hvib_re)
    assert expected_result == round(computed_result.real, 5) + round(computed_result.imag,5) * 1j





