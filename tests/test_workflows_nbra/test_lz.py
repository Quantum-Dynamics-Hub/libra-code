"""
Unit and regression test for the workflows/nbra module in the Libra package
"""

from libra_py.workflows.nbra import lz
from libra_py import units
import pytest
import math
import sys
import os
 
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *




def generate_hvib_data(model):
    """
    This function generates model data to be used to test the functions in the libra_py.workflows.nbra.lz module
    """
    dt     = 0.01  # atomic units
    nsteps = 1000
    hvib   = []    
    if model == 1:

        # This model consists of one state with constant zero energy and another which is a sin function. Therefore,
        # the energy gap is the same sin funciton as well. The model is designed such that the minima of the sin function
        # occurs at poitns which have energy gaps of 0.1 eV. Such minima occur on at t = 2*math.pi*n + 3*math.pi*0.5, 
        # for integer values of n
 
        de     = 0.1*units.ev2Ha
        for t in range(nsteps):
            hvib.append(CMATRIX(2,2))
            hvib[t].set(0,0, -0.0+0.0j);        hvib[t].set(0,1, 0.0*(0.0+1.0j));
            hvib[t].set(1,0,  0.0*(0.0-1.0j));  hvib[t].set(1,1, 2*de + de*math.sin(t*dt)+0.0j);


    return hvib




def generate_params(model):
    """
    This function generates the parameters dictionary to be used to test the functions in the libra_py.workflows.nbra.lz module
    """
    params = {}
    if model == 1:
        params["dt"]                 = 0.01
        params["T"]                  = 300.0
        params["Boltz_opt_BL"]       = 1  
        params["gap_min_exception"]  = 0
        params["target_space"]       = 1

    if model == "test_adjust_SD_probabilities_1":
        params["excitations"] = [ [ [ [1, 2, -5, -6], [1, 3, -5, -6], [1, 4, -5, -6], [3, 2, -5, -6], [4, 2, -5, -6] ] ] ]

    #if model == "test_adjust_SD_probabilities_2":
    #    params["excitations"] = [ [ [ [1, 2], [1, 3], [1, 4], [3, 2], [4, 2] ] ] ]

    return params



@pytest.mark.parametrize("hvib, params, expected_result", [
    ( generate_hvib_data(1), generate_params(1), 0.9942 ),   
])
def test_Belyaev_Lebedev(hvib, params, expected_result):
    """Tests that the LZ probabilities are computed correctly according to the NBRA BLSH scheme"""   
    P = lz.Belyaev_Lebedev(hvib, params)

    computed_result = 0
    for t in range(len(P)):
        if P[t].get(0,0) != 1:
            computed_result = P[t].get(0,1)
            assert expected_result == round(computed_result, 4)
            break



@pytest.mark.parametrize("Hvib, params, expected_result", [
    ( [generate_hvib_data(1)], generate_params(1), 0.9942 ),
])
def test_lz(Hvib, params, expected_result):
    """Tests that the LZ probabilities are computed correctly according to the NBRA BLSH scheme"""
    params["istate"]             = 1                  # From 0
    params["T"]                  = 300.0              # Temperature, K
    params["target_space"]       = 1
    params["gap_min_exception"]  = 0
    params["Boltz_opt_BL"]       = 1                  # Option to incorporate hte frustrated hops into BL probabilities
    params["outfile"]            = "_out_Markov_.txt" # output file
    params["evolve_Markov"]      = True               # Rely on the Markov approach
    params["evolve_TSH"]         = False              # don't care about TSH
    params["ntraj"]              = 1                  # how many stochastic trajectories
    params["init_times"]         = [0]                # starting points for sub-trajectories
    params["do_output"]          = True               # request to print the results into a file
    params["do_return"]          = False              # request to not store the date in the operating memory    
    params["return_probabilities"] = True

    res, P = lz.run(Hvib, params)

    computed_result = 0
    for t in range(len(P[0])):
        if P[0][t].get(0,0) != 1:
            computed_result = P[0][t].get(0,1); #print (computed_result); sys.exit(0)
            assert expected_result == round(computed_result, 4)
            break
   



def get_data(filename, data_dim):
    X = MATRIX(data_dim, data_dim)
    X.Load_Matrix_From_File(filename)
    return X




@pytest.mark.parametrize("P, params, expected_result", [
    ( [[get_data("supporting_files/5_state_LZ_probabilities.txt", 5)]], generate_params("test_adjust_SD_probabilities_1"), [[get_data("supporting_files/5_state_LZ_probabilities_expected.txt", 5)]] ),
])
def test_adjust_SD_probabilities(P, params, expected_result):
    """Tests that the LZ probabilities are computed correctly according to the NBRA BLSH scheme"""
    computed_result = lz.adjust_SD_probabilities(P, params)  
    assert expected_result[0][0].get(1,4) == computed_result[0][0].get(4,1)




