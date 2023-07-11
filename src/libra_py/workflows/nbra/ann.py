#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
"""

.. module:: ann
   :platform: Unix, Windows
   :synopsis: This module implements routines to train and use ANN-based Hamiltonian models for NBRA dynamics
.. moduleauthor:: Alexey V. Akimov


"""

import sys
import cmath
import math
import os
import multiprocessing as mp
#import multiprocess as mp
import time
import json
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import libra_py.common_utils as comn
import util.libutil as comn

import libra_py.data_read as data_read
import libra_py.data_conv as data_conv
from . import decoherence_times as dectim
#import libra_py.tsh as tsh
import libra_py.tsh_stat as tsh_stat
import libra_py.units as units
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.fit as fit
import libra_py.influence_spectrum as infsp

class tmp:
    pass

def preprocess_data(npdata, verbose=False):
    """
    Args:
        data  ( numpy array of floats ): original data
        
    Return:
        nparray, float, float: transformed data, scaling factor, shift factor
        
        transformed_data = scaling_factor * (original_data - shift_factor)
    """
    
    # convert to the np array
    #npdata = np.array( data )
    
    # Compute average
    data_ave = np.average(npdata)
    if verbose:
        print(F"Average of the original data set = {data_ave}")

    # Center the data
    npdata = npdata - data_ave
    
    if verbose:
        tmp = np.average(npdata)
        print(F"Average of the shifted data set = {tmp}")

    min_val = np.amin(npdata)
    max_val = np.amax(npdata)
    if verbose:
        print(F"minimum = {min_val}  maximum = {max_val}")

    scaling = (0.75 / (max_val - min_val) ) 
    npdata =  scaling * npdata
    
    return npdata, scaling, data_ave
    

def find_maxima(w, data, fraction_tol):
    """
    A generic function, in principle - it looks for the maxima 
    in the given data. It returns all the maxima and a certain 
    fraction of the leading maxima (highest values) as determined
    by the threshold parameter

    Args:
        w ( list of floats ): x axis values, e.g. frequencies
        data ( list of floats ): y axis values, e.g. intencities
        fraction_tol ( float ): the fraction of the leading frequencies to select

    Returns:
        list of floats, list of floats: the first list are all maxima, the second - only the leading ones
    """

    res = []
    n = len(data)
    
    for i in range(1,n-1):
        if data[i] > data[i+1] and data[i] > data[i-1]:
            res.append([ i, data[i]])
            
    res_srtd = merge_sort(res)
    
    n = len(res_srtd)
    
    final_res = []
    for i in range(n):
        I = res_srtd[n-1-i][0]
        j2 = res_srtd[n-1-i][1]
        final_res.append( [w[I], j2] )
        
    
    max_j2 = final_res[0][1]
    
    leading_freqs = [ final_res[0][0] ]
    for i in range(1,n):
        if (final_res[i][1]/ final_res[0][1]) > fraction_tol:
            leading_freqs.append( final_res[i][0] )
                
    return final_res, leading_freqs



def make_input(t, w, tau):
    """
    Inputs to the ANN-based Hamiltonians for NBRA:

    This function converts time into a set of inputs (periodic and exponenetial functions of time)
    to the ANN
    
    The input to the ANN is a (dim1+dim2)-dimensional vector: 

    (sin(w_0*t), sin(w_1*t), ... sin(w_{dim1-1}*t), exp(-t/tau_0), exp(-t/tau_1), ..., exp(-t/tau_{dim2-1}) )

    Args:
        t ( float ): time [units: a.u. of time]
        w ( list of floats ): frequencies of the ANN modes [units: cm^-1]
        tau ( list of floats ): the exponential decay parameters  

    """
    
    dim1 = len(w)
    dim2 = len(tau)
    
    inputs  = MATRIX(dim1+dim2, 1)    
            
    # Inputs
    for idim in range(dim1):
        inputs.set(idim, 0, math.sin(w[idim] * units.wavn2au * t) )
    for idim in range(dim2):
        inputs.set(dim1+idim, 0, math.exp(-t/tau[idim]) )
            
    return inputs



def make_training_data(_X, indices_subset, _w_lead, tau, dt=1.0 * units.fs2au, deriv_scaling=50.0):
    """
    Args:        
        _X ( np.array of floats ): the whole dataset - energy gaps or couplings
        indices_subset ( list of ints ): the indices of the datapoints from the entire dataset to use 
        w_lead ( list of floats ): dominant frequencies in the data set [ units: cm^-1 ]
        tau ( list of floats ): exponential decay parameters for the ANN inputs [ units: a.u.]
        dt ( float ): the timestep between the data points [ units: a.u. ]
        deriv_scaling  ( float ): parameter that weights the finite differences in the error functions 
        
    Returns:
        MATRIX(dim1+dim2, n_samples), MATRIX(2, n_samples): the inputs, the targets
        
    """
    
    n_samples = len(indices_subset)
    
    dim1 = len(_w_lead)
    dim2 = len(tau)

    print(F"dimentionality of input is {dim1 + dim2}")

    inputs  = MATRIX(dim1 + dim2, n_samples)
    targets = MATRIX(2, n_samples)   # values of the function and the "scaled derivatives"
    
    for istep, step in enumerate(indices_subset):
        t = dt * step  # in a.u. now
        
        # Inputs
        for idim in range(dim1):
            inputs.set(idim, istep, math.sin(_w_lead[idim] * units.wavn2au * t) )
        for idim in range(dim2):
            inputs.set(dim1+idim, istep, math.exp(-t/tau[idim]) )

            
        # Targets
        ### Function
        targets.set(0, istep, _X[step])
        
        ### Scaled derivative 
        if step>0 and step<len(indices_subset)-1:
            targets.set(1, istep, 0.5*deriv_scaling*(_X[step+1] - _X[step-1]) )
        elif step==0:
            targets.set(1, istep, deriv_scaling*(_X[step+1] - _X[step] ) )
        elif step==len(indices_subset)-1:
            targets.set(1, istep, deriv_scaling*(_X[step] - _X[step-1] ) )

        
    return dim1+dim2, inputs, targets




def step1_workflow(data,  _params, _plt):
    """
    Args:
        * data ( numpy array of floats ): original data
        * _params ( dict )
            - dt ( float ): FT timestep [ units: fs, default: 1.0 ]
            - wspan ( float ): range of frequencies to expect [ units: cm^-1, default : 500.0 ]
            - dw ( float ): frequency scale resolution [ units: cm^-1, default : 1.0 ]
            - do_output ( Bool ): whether to printout extra info [ default : False ]
            - do_center ( Bool ): whether to center the data by cubtracting its average [ default : True ]
            - acf_type ( int ): selection of the ACF type [ default : 1]
            - data_type ( int ): selection of the data type [ defualt : 0 ]
            - leading_w_fraction ( float ): a fraction of the leading frequencies to be used in the ANN training [ default: 0.0001]
            - deriv_scaling ( float ): the weight of the slopes in the final error function [ default: 50.0 ]
            - tau ( list of floats ): the exponential decay parameters for the ANN inputs [ units: a.u. of time, default: [5000.0] ]
            - output_files_prefix ( string ): a common prefix for all filex generated by this function
        * plt ( gnuplot plotting instance )
    """

    params = dict(_params)

    critical_params = []
    default_params = { "dt": 1.0, "wspan":500.0, "dw":1.0, "do_output":False, "do_center":True, 
                       "acf_type":1, "data_type":0, 
                       "leading_w_fraction":0.0001,  "deriv_scaling":50.0,
                       "tau":[5000.0], "training_steps":[-1],
                       "output_files_prefix":"data"
                     }
    comn.check_input(params, default_params, critical_params)


    w_fract = params["leading_w_fraction"]
    dt = params["dt"] * units.fs2au
    deriv_scaling = params["deriv_scaling"]
    output_files_prefix = params["output_files_prefix"]
    tau = params["tau"]
    training_steps = params["training_steps"]
    if training_steps == [-1]:
        training_steps = list(range(len(data)))
    


    ## ======== Data shaping ============

    data_transformed, data_scaling, data_shift = preprocess_data(data, True) 


    ## ======== Forier Transform =======
    all_data = []

    for idat, dat in enumerate(data_transformed):
        x = MATRIX(1,1)    
        x.set(0,0, dat)
        all_data.append(x)
        
    T, ACF, uACF, W, J, J2 = infsp.recipe1(all_data, params)


    ## ============ Plot the spectra =============

    figure = _plt.figure(num=None, figsize=(3.21*2, 2.41), dpi=300, edgecolor='black', frameon=True)
 
    _plt.subplot(1,2,1)
    _plt.xlabel('Time, fs',fontsize=10)
    _plt.ylabel('normalized ACF',fontsize=10)
    _plt.plot(T, ACF, color="blue", label="", linewidth=2)

    _plt.subplot(1,2,2)
    _plt.xlabel('Wavenumber, $cm^{-1}$',fontsize=10)
    _plt.ylabel('Influence spectrum',fontsize=10)
    _plt.plot(W, J2, color="blue", label=F"", linewidth=2)

    _plt.tight_layout()
    _plt.savefig(F"{output_files_prefix}-acf-ifs.png", dpi=300)
    _plt.show()


    ## ============= Find the leading terms ===============
    z, w_lead = find_maxima(W, J2, w_fract)  # energy gap


 
    ## ============ Make training data ====================
    ann_input_dim, training_input, training_target = make_training_data(data_transformed, training_steps, w_lead, tau, dt, deriv_scaling )


    ## ============= Now save meta-information and the data ==========
    meta_info = { "ann_input_dim": ann_input_dim, "w_lead":w_lead, "tau":tau, 
                  "data_scaling":data_scaling, "data_shift":data_shift,
                  "training_steps":training_steps,
                  "training_input_n_rows":training_input.num_of_rows, "training_input_n_cols":training_input.num_of_cols,
                  "training_target_n_rows":training_target.num_of_rows, "training_target_n_cols":training_target.num_of_cols
                }

    with open(F"{output_files_prefix}-meta.json", 'w') as outfile:
        json.dump(meta_info, outfile)

    
    training_input.bin_dump(F"{output_files_prefix}-training-input.MATRIX")
    training_target.bin_dump(F"{output_files_prefix}-training-target.MATRIX")




def plot_error(_plt, err, filename="ANN_err.png"  ):
    """
    Args: 
        _plt (matplotlib instance)
        err (doubleList ): error values for all iterations

    """
    
    y = Cpp2Py(err)
    x = list( range(len(y)))
       
    _plt.rc('axes', titlesize=12)      # fontsize of the axes title
    _plt.rc('axes', labelsize=12)      # fontsize of the x and y labels
    _plt.rc('legend', fontsize=10)     # legend fontsize
    _plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
    _plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
    _plt.rc('figure.subplot', left=0.2)
    _plt.rc('figure.subplot', right=0.95)
    _plt.rc('figure.subplot', bottom=0.13)
    _plt.rc('figure.subplot', top=0.88)
    
    
    figure = _plt.figure(num=None, figsize=(3.21*1, 2.41), dpi=300, edgecolor='black', frameon=True)

    _plt.subplot(1,1,1)
    _plt.xlabel('Iteration number',fontsize=10)
    _plt.ylabel('Error',fontsize=10)
    _plt.plot(x, y, color="black", label="", linewidth=2)
    _plt.legend(fontsize=6.75, ncol=1, loc='upper left')
        
    _plt.tight_layout()
    _plt.savefig(filename, dpi=300)
    _plt.show()
    
    

def plot_training(_plt, _ann, _training_input, _training_target, mode1=0, mode2=2, 
                  filename="ANN_training.png"  ):
    """
    Args: 
        _plt (matplotlib instance)
        _ann (NeuralNetwork object)
        _training_input ( MATRIX(n_inputs, n_patterns )
        _training_target ( MATRIX(n_outputs, n_patterns )

    """
    
    last_indx = _ann.Nlayers - 1 
        
    
    # Predict
    # There is only 1 output in this case - we take the first index
    training_out = _ann.propagate(_training_input)
    y_predict = data_conv.MATRIX2nparray(training_out[last_indx])[0, :]
    
    #print(F"len(training_out) = {len(training_out)}")
    #print(F"len(y_predict) = {len(y_predict)}")

    # Given targets
    # Targets also have only 1 components in this case
    y_targets = data_conv.MATRIX2nparray(_training_target)[0, :]
    
    #print(F"len(y_targets) = {len(y_targets)}")

    nsteps = len(y_predict)    # one less, so it would work for NACs too
    steps_train = list( range(nsteps) ) 
    
    #print(F"nsteps = {nsteps}")
    #print(F"steps_trian = {steps_train}")

    
    x_mode1 = data_conv.MATRIX2nparray(_training_input)[mode1, :]
    
    x_mode2 = data_conv.MATRIX2nparray(_training_input)[mode2, :]
    
    #print(x_mode1)
    #print(x_mode2)
    
    _plt.rc('axes', titlesize=12)      # fontsize of the axes title
    _plt.rc('axes', labelsize=12)      # fontsize of the x and y labels
    _plt.rc('legend', fontsize=10)     # legend fontsize
    _plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
    _plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
    _plt.rc('figure.subplot', left=0.2)
    _plt.rc('figure.subplot', right=0.95)
    _plt.rc('figure.subplot', bottom=0.13)
    _plt.rc('figure.subplot', top=0.88)
    
    
    figure = _plt.figure(num=None, figsize=(3.21*3, 2.41), dpi=300, edgecolor='black', frameon=True)

    _plt.subplot(1,3,1)
    _plt.xlabel('Timestep',fontsize=10)
    _plt.ylabel('Property',fontsize=10)
    _plt.plot(steps_train, y_predict[:nsteps], color="black", label="ANN", linewidth=2)
    _plt.plot(steps_train, y_targets[:nsteps], color="green", label="reference", linewidth=2)
    _plt.legend(fontsize=6.75, ncol=1, loc='upper left')
    

    _plt.subplot(1,3,2)
    _plt.xlabel(F'Mode {mode1}',fontsize=10)
    _plt.ylabel('Property',fontsize=10)    
    x_mode1 = data_conv.MATRIX2nparray(_training_input)[mode1, :]
    _plt.plot(x_mode1[:nsteps], y_predict[:nsteps], color="black", label="ANN", linewidth=2)
    _plt.plot(x_mode1[:nsteps], y_targets[:nsteps], color="green", label="reference", linewidth=2)
    _plt.legend(fontsize=6.75, ncol=1, loc='upper left')
    
    _plt.subplot(1,3,3)
    _plt.xlabel(F'Mode {mode2}',fontsize=10)
    _plt.ylabel('Property',fontsize=10)    
    x_mode2 = data_conv.MATRIX2nparray(_training_input)[mode2, :]
    _plt.plot(x_mode2[:nsteps], y_predict[:nsteps], color="black", label="ANN", linewidth=2)
    _plt.plot(x_mode2[:nsteps], y_targets[:nsteps], color="green", label="reference", linewidth=2)
    _plt.legend(fontsize=6.75, ncol=1, loc='upper left')

    
    _plt.tight_layout()
    _plt.savefig(filename, dpi=300)
    _plt.show()
    
    

def step2_workflow(_params, _ann_params, rnd, _plt):
    """
    Loads the meta-info about the data + training and target data from the files and trains ANN

    Args:
        _params ( dict ): general control parameters
            - input_files_prefix ( string ): the common prefix for various files related to the same data set [ default: "data"]
            - ann_arch ( list of ints ): architecture of the hidden layers [ default: [10, 10] ]
        _ann_params ( dict ): parameters controlling the approach to train the ANNs
        rnd ( Random ): random number generator object
        _plt ( matplotlib instance )

    """

    params = dict(_params)

    critical_params = []
    default_params = { "input_files_prefix":"data", "ann_arch":[10,10] , "mode1":0, "mode2":1   }
    comn.check_input(params, default_params, critical_params)
   
    input_files_prefix = params["input_files_prefix"]
    ann_arch = params["ann_arch"]
    mode1 = params["mode1"]
    mode2 = params["mode2"]


    ## ================= Load the files ===========
    json_file = open(F"{input_files_prefix}-meta.json")
    meta_info = dict(json.load(json_file))
    json_file.close()

    n,m = meta_info["training_input_n_rows"], meta_info["training_input_n_cols"]
    training_input = MATRIX(n,m)
    training_input.bin_load(F"{input_files_prefix}-training-input.MATRIX")

    n,m = meta_info["training_target_n_rows"], meta_info["training_target_n_cols"]
    training_target = MATRIX(n,m)
    training_target.bin_load(F"{input_files_prefix}-training-target.MATRIX")


    ## ======================= Train the ANN ==================
    input_dim = meta_info["ann_input_dim"]
 
    
    ANN = NeuralNetwork( Py2Cpp_int( [input_dim] + ann_arch + [2] ) ) 
    ANN.init_weights_biases_uniform(rnd, -1.1, 1.1, -1.1, 1.1)


    ann_params = dict(_ann_params)
    ann_critical_params = []
    ann_default_params = { "num_epochs":25, 
                       "steps_per_epoch":1000, 
                       "error_collect_frequency":1,
                       "epoch_size":1, 
                       "learning_method":2,
                       "learning_rate":0.001,
                       "weight_decay_lambda":0.000,
                       "etha":1.0,
                       "momentum_term": 0.75,
                       "verbosity":1,
                       "a_plus":1.1, "a_minus":0.75,
                       "dB_min":0.000001, "dB_max":10.1,
                       "dW_min":0.000001, "dW_max":10.1,
                     }
    comn.check_input(ann_params, ann_default_params, ann_critical_params)

    err = ANN.train(rnd, ann_params, training_input, training_target )
    ANN.save(F"{input_files_prefix}-ann.xml")


    ## =============== Plotting ==================

    plot_error(_plt, err, F"{input_files_prefix}-ann-error.png")

    plot_training(_plt, ANN, training_input, training_target, mode1, mode2, F"{input_files_prefix}-ann-training.png")


def load_ann_and_parameters(params, nstates=2, prefix="./"):
    """
    This is a convenience function that loads the ANNs and meta info for all the data sets

    Args:
        params ( dict ): this dictionary is to be updated with what is read from the files
        nstates ( int ): the number of states for which the ANNs are ready - the number and the names 
            of the files that are present should be consistent with this number [ default: 2]
        prefix ( string ): the location of all the files [ default: "./"]

    Returns:
        None: but updates the `params` dictionary that is provided as the input

    """

    ## =================== ENERGIES ====================
    params.update( { "gap_ann":[] , "gap_ann_w":[], "gap_ann_tau":[], "gap_ann_shift":[], "gap_ann_scale":[] } )

    for i in range(1, nstates):

        params["gap_ann"].append([]) 
        params["gap_ann_w"].append([]) 
        params["gap_ann_tau"].append([]) 
        params["gap_ann_shift"].append([]) 
        params["gap_ann_scale"].append([]) 

        ## ================= Load the meta-info ===========
        json_file = open(F"{prefix}e{i-1}{i}-meta.json")
        meta_info = dict(json.load(json_file))
        json_file.close()

        ## ================= Load the data ================
        params["gap_ann"][i-1] = NeuralNetwork( F"{prefix}e{i-1}{i}-ann.xml")
        params["gap_ann_w"][i-1] =  list(meta_info["w_lead"]) 
        params["gap_ann_tau"][i-1] = list(meta_info["tau"]) 
        params["gap_ann_shift"][i-1] = meta_info["data_shift"]
        params["gap_ann_scale"][i-1] = 1.0/meta_info["data_scaling"]

        

    ## =================== COUPLINGS ====================
    params.update( { "nac_ann":[] , "nac_ann_w":[], "nac_ann_tau":[], "nac_ann_shift":[], "nac_ann_scale":[] } )


    cnt = 0
    for i in range(0, nstates):
        for j in range(i+1, nstates):

            params["nac_ann"].append([]) 
            params["nac_ann_w"].append([]) 
            params["nac_ann_tau"].append([]) 
            params["nac_ann_shift"].append([]) 
            params["nac_ann_scale"].append([]) 


            ## ================= Load the meta-info ===========
            json_file = open(F"{prefix}nac{i}{j}-meta.json")
            meta_info = dict(json.load(json_file))
            json_file.close()

            ## ================= Load the data ================
            params["nac_ann"][cnt] =  NeuralNetwork( F"{prefix}nac{i}{j}-ann.xml") 
            params["nac_ann_w"][cnt] =  list(meta_info["w_lead"]) 
            params["nac_ann_tau"][cnt] = list(meta_info["tau"]) 
            params["nac_ann_shift"][cnt] = meta_info["data_shift"]
            params["nac_ann_scale"][cnt] = 1.0/meta_info["data_scaling"]

            cnt = cnt + 1




def compute_model_nbra_ann(q, params, full_id):
    """   
    This function construct the ANN-based Hamiltonians for NBRA calculations

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof, but they do not really affect anything
        params ( dictionary ): model parameters

            * **params["timestep"]** ( int ):  index of time step of this trajectory 
                this parameter is varied by the dynamics function that calls this one
            * **params["istep"]** ( int ):  index of time step of the initial condition, together with
                the `timestep` parameter, it determines the global time since the time origin 
            * **params["dt"]** ( dt ): magnitude of the time step [units: a.u. of time] 
            * **params["nstates"]** ( int ): the total number of states

            * **params["gap_ann"]** ( list of nstates-1 NeuralNetwork objects ): ANNs for E_i - E_{i-1} gaps, i = 1, ... nstates
            * **params["gap_ann_w"]** ( list of nstates-1 lists of floats ): frequencies needed to 
                construct the inputs to the ANN for  E_i - E_{i-1} gaps, i = 1, ... nstates, the number of frequencies can vary 
                depending on the gap we describe
            * **params["gap_ann_tau"]** ( list of nstates-1 lists of floats ): exponential decay parameter needed to 
                construct the inputs to the ANN for  E_i - E_{i-1} gaps, i = 1, ... nstates, the number of frequencies can vary 
                depending on the gap we describe
            * **params["gap_ann_scale"]** ( list of nstates-1 floats): scaling parameters for the output produced by the ANN to convert 
                it to the desired gap
            * **params["gap_ann_shift"]** ( list of nstates-1 floats): shift parameters for the output produced by the ANN to convert 
                it to the desired gap

            * **params["nac_ann"]** ( list of N(N-1)/2 NeuralNetwork objects ): ANNs for NAC_ij, j>i gaps, i = 1, ... nstates (N = nstates) 
            * **params["nac_ann_w"]** ( list of N(N-1)/2 lists of floats ): frequencies needed to 
                construct the inputs to the ANN for NAC_ij, j>i gaps, i = 1, ... nstates (N = nstates), the number of frequencies can vary 
                depending on the nac we describe
            * **params["nac_ann_tau"]** ( list of N(N-1)/2 lists of floats ): exponential decay parameter needed to 
                construct the inputs to the ANN for NAC_ij, j>i gaps, i = 1, ... nstates (N = nstates), the number of frequencies can vary 
                depending on the nac we describe
            * **params["nac_ann_scale"]** ( list of N(N-1)/2 floats): scaling parameters for the output produced by the ANN to convert 
                it to the desired nac
            * **params["nac_ann_shift"]** ( list of N(N-1)/2 floats): shift parameters for the output produced by the ANN to convert 
                it to the desired nac

        
    Returns:       
        PyObject: obj, with the members:

            * obj.ham_adi ( CMATRIX(n,n) ): adiabatic electronic Hamiltonian 
            * obj.nac_adi ( CMATRIX(n,n) ): scalar non-adiabatic couplings
            * obj.hvib_adi ( CMATRIX(n,n) ): adiabatic vibronic Hamiltonian 
            * obj.basis_transform ( CMATRIX(n,n) ): adi-to-dia transformation matrix, identity here
            * obj.time_overlap_adi ( CMATRIX(n,n) ): time-overlap matrix, identity here

            
    """

    hvib_adi, basis_transform, time_overlap_adi = None, None, None
  
    Id = Cpp2Py(full_id)
    indx = Id[-1]    
    timestep = params["timestep"]
    istep = params["istep"]
    dt = params["dt"]
    nadi = params["nstates"]
    
    #print(F"timestep = {timestep} istep = {istep}  dt = {dt} nadi = {nadi}")
        
    t = dt * (timestep + istep)
            
    #============ Electronic Hamiltonian =========== 
    ham_adi = CMATRIX(nadi, nadi)         
    
    
    ei = 0.0 
    for i in range(1, nadi):
        ann = params["gap_ann"][i-1]          # ANN for gap between states i and i-1
        w = params["gap_ann_w"][i-1]          # frequencies
        tau = params["gap_ann_tau"][i-1]      # tau parameters
                
        ann_input = make_input(t, w, tau)
        last_indx = ann.Nlayers - 1

        de = ann.propagate(ann_input)[last_indx].get(0,0)   
        
        gap_scale = params["gap_ann_scale"][i-1]
        gap_shift = params["gap_ann_shift"][i-1]
        
        #print(F"w = {w} tau = {tau} de = {de}  gap_scale = {gap_scale}  gap_shift = {gap_shift}")
        
        de = gap_scale * de + gap_shift
        
        ei = ei + de
        ham_adi.set(i, i,  ei * (1.0+0.0j) ) 
    
        
    
    #============ NAC ===========        
    nac_adi = CMATRIX(nadi, nadi)    
        
    cnt = 0
    for i in range(0, nadi):
        for j in range(i+1, nadi):
            ann = params["nac_ann"][cnt]     # ANN for NAC between states i and j
            w = params["nac_ann_w"][cnt]     # frequencies
            tau = params["nac_ann_tau"][cnt] # tau parameter
                
            ann_input = make_input(t, w, tau)
            last_indx = ann.Nlayers - 1
        
            nac = ann.propagate(ann_input)[last_indx].get(0,0)
        
            nac_scale = params["nac_ann_scale"][cnt]
            nac_shift = params["nac_ann_shift"][cnt]
        
            nac = nac_scale * nac + nac_shift

            #print(F"w = {w} tau = {tau} nac = {nac}  nac_scale = {nac_scale}  nac_shift = {gap_shift}")
                
            nac_adi.set(i,j,  nac * (1.0+0.0j) ) 
            nac_adi.set(j,i,  nac * (-1.0+0.0j) ) 
        
            cnt = cnt + 1
          
    
    #============ Vibronic Hamiltonian ===========        
    hvib_adi = CMATRIX(nadi, nadi)
    hvib_adi = ham_adi - (0.0+1j) * nac_adi
        
    #=========== Basis transform, if available =====
    basis_transform = CMATRIX(nadi, nadi)            
    basis_transform.identity()        
                                                
    #========= Time-overlap matrices ===================
    time_overlap_adi = CMATRIX(nadi, nadi)            
    time_overlap_adi.identity()    
    
    
    obj = tmp()
    obj.ham_adi = ham_adi
    obj.nac_adi = nac_adi
    obj.hvib_adi = hvib_adi
    obj.basis_transform = basis_transform
    obj.time_overlap_adi = time_overlap_adi
            
    return obj


class tmp1:
    pass

def ann_getstate(self):    
    x = tmp1()
    
    x.B = self.B
    x.grad_b = self.grad_b
    x.grad_b_old = self.grad_b_old
    x.dB = self.dB
    x.dBold = self.dBold
    x.W = self.W
    x.grad_w = self.grad_w
    x.grad_w_old = self.grad_w_old
    x.dW = self.dW
    x.dWold = self.dWold
    x.Nlayers = self.Nlayers 
    x.Npe = self.Npe
    x.sz_x = self.sz_x
    x.sz_y = self.sz_y
    
    return (x, )

def ann_setstate(self, x):
    
    self.allocate(x[0].Npe)

    self.B = x[0].B
    self.grad_b = x[0].grad_b
    self.grad_b_old = x[0].grad_b_old
    self.dB = x[0].dB
    self.dBold = x[0].dBold
    self.W = x[0].W
    self.grad_w = x[0].grad_w
    self.grad_w_old = x[0].grad_w_old
    self.dW = x[0].dW
    self.dWold = x[0].dWold
    self.Nlayers = x[0].Nlayers 
    self.Npe = x[0].Npe
    self.sz_x = x[0].sz_x
    self.sz_y = x[0].sz_y

#NeuralNetwork.__getinitargs__ = ann_getinitargs
NeuralNetwork.__getstate__ = ann_getstate
NeuralNetwork.__setstate__ = ann_setstate


