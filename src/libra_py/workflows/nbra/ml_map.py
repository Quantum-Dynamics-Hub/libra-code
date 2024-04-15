#***********************************************************
# * Copyright (C) 2024 Mohammad Shakiba and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

"""
.. module:: ml_map
   :platform: Unix
   :synopsis: This module creates models to map KS Hamiltonian matrices from
              one level of theory to another and compute the properties of the
              ML generated KS Hamiltonian.

.. moduleauthor:: Mohammad Shakiba and Alexey V. Akimov

"""

import os
import sys
import glob
import time
import numpy as np
import multiprocessing as mp
import scipy.sparse as sp
import pickle
import joblib
# For future works!
# import tensorflow as tf
# from tensorflow.keras import layers, models
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_squared_error, accuracy_score, mean_absolute_error, r2_score
from liblibra_core import *
import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.packages.dftbplus.methods as DFTB_methods
from libra_py import units, molden_methods, data_conv
import util.libutil as comn


def find_indices(params):
    """
    This function finds the indices of a set of file in a directory
    sorts them and return them in a list
    Args:
        params (dict): The details of the parameters in the dictionary are given in the
                       train function.
    Returns:
        indices (list): The list of indices of the reference calculation
    """
    # The indices are generated based on the output matrices
    files = glob.glob(f'{params["path_to_output_mats"]}/*{params["output_property"]}*npy')
    indices = []
    for file in files:
        # Generate the index in a file
        index = file.split('/')[-1].split('_')[-1].replace('.npy','')
        #print(index)
        indices.append(int(index))

    # Now sort the indices since they will be used
    indices = np.sort(indices)
    return indices


def read_data(params, path, prefix):
    """
    This function is used to read the input and output (reference) matrices
    from the path it is given to in the params variable
    It is not a general function but rather an auxiliary one so it's not
    suitable in data_read module in libra_py.
    Args:
        params (dict): The details of the parameters in the dictionary are given in the
                       train function.
        path (string): The path to data
        prefix (string): The prefix used to save the data
    Returns:
        data (nparray): The read data in the path
    """
    data = []
    for i in params["train_indices"]:
        #tmp = np.load(F'{path}/{prefix}_{params["input_property"]}_{i}.npy')
        tmp = np.load(F'{path}/{prefix}_{i}.npy')
        data.append(tmp)
    data = np.array(data)
    return data 

def partition_matrix(params, matrix):
    """
    This function is used to partition an input matrix
    There are different partitioning approaches proposed in [Ref]
    Here, we implement these methods:
    1- 'eqaul': Equal splitting of the vectorized upper triangular part of the matrix
    2- 'block': We find the blocks related to interaction of two atoms in the matrix
                and then vectorize it. We first do this for off-diagonal blocks and
                then for the diagonal blocks (interaction of each atom with itself)
    3- 'atomic': The section of a matrix which shows the interaction of each 
                   atom with all other atoms
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        matrix (nparray): The matrix to be partitioned
    Returns:
        partitioned_matrix (list): A list containing the partitions of the matrix
    """
    # =========== 1- 'equal' partitioning
    if params["partitioning_method"]=="equal":
        # generate the upper indices of the matrix
        upper_indices = np.triu_indices(matrix.shape[0])
        # Extract the upper matrix as a vector
        upper_vector = matrix[upper_indices]
        npartition = params["npartition"]
        partition_points = np.linspace(0, upper_vector.shape[0], npartition+1, dtype=int, endpoint=True)
        #print(partition_points, upper_vector.shape[0])
        partitioned_matrix = []
        for i in range(len(partition_points)-1):
            start = partition_points[i]
            #if i==len(partition_points)
            end = partition_points[i+1]
            partitioned_matrix.append(upper_vector[start:end])
    # =========== 2- 'block' partitioning
    elif params["partitioning_method"]=="block":
        print("To be implemented")
    # =========== 3- 'atomic' partitioning
    elif params["partitioning_method"]=="atomic":
        #print("To be implemented")
        # Find the Log files where the AO matrices are stored
        # From these files, one can recognize the indices 
        # corresponding to each atom in the matrix 
        AO_Log_files = glob.glob(f"{params['path_to_sample_files']}/{params['prefix']}*.Log")
        if params["input_partition"]:
            for i in range(len(AO_Log_files)):
                if not "_ref" in AO_Log_files[i]:
                    ao_file = AO_Log_files[i]
                    break
        else:
             for i in range(len(AO_Log_files)):
                if "_ref" in AO_Log_files[i]:
                    ao_file = AO_Log_files[i]
                    break
        ao_indices = CP2K_methods.atom_components_cp2k(ao_file)
        mat_size = matrix.shape[0]
        partitioned_matrix = []
        for i in range(len(ao_indices)):
            partition = []
            for indx1 in ao_indices[i]:
                for indx2 in range(indx1, mat_size):
                    partition.append(matrix[indx1, indx2])
            partitioned_matrix.append(partition)

    return partitioned_matrix

def partition_data(params, data):
    """
    This function uses the partition_matrix to partition a set of data
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        data (list): A list of numpy arrays of matrices
    Returns:
        partitioned_data (list): A list containing the partitioned matrices
    """
    partitioned_data = []
    #print(data.shape)
    for i in range(len(data)):
        tmp = partition_matrix(params, data[i])
        partitioned_data.append(tmp)
    #partitioned_data = np.array(partitioned_data)
    #print(len(partitioned_data), len(partitioned_data[0]), len(partitioned_data[0][0]))
    return partitioned_data

def scale_partition(params, partition):
    """
    This function is used to scale the 'set' of a specific partition of
    the matrix based on the scalers in scikit-learn package
    It returns both the scaler and the scaled data
    The scaler is needed for transforming back the data
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        partition (list): The partition of a matrix
    Returns:
        scaler (scikit-learn scaler): The scaler function of scikit-learn
        scaled_data (nparra): The scaled partition
    """
    if params["scaler"].lower()=="standard_scaler":
        scaler = StandardScaler()     
    elif params["scaler"].lower()=="minmax":
        scaler = MinMaxScaler()
    scaler.fit(partition)
    scaled_data = scaler.transform(partition)
    #print(partition.shape)
    return (scaler, scaled_data)

def scale_data(params, data):
    """
    This function is used to scale all the partitions in the inputs or outputs
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        scaled_data (list): The list of data to be partitioned
    Returns:
        scalers (list): A list of scalers
        scaled_data (list): A list of scaled data 
    """
    scalers = []
    scaled_data = []
    for i in range(len(data[0])): # The number of partitioned and available vectors
        p1 = []
        for k in range(len(data)):
            p1.append(data[k][i])
        p1 = np.array(p1)
        #print('In for loop of scale:', p1.shape)
        res = scale_partition(params, p1)
        scalers.append(res[0])
        scaled_data.append(res[1])
    #return scalers, np.array(scaled_data)
    return (scalers, scaled_data)

def train_partition(params, input_data, output_data, output_scaler):
    """
    This function builds and trains a KRR model
    using input_data and output_data
    We will use the scaled input and output data to train the model
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        input_data (list): The input data vectors for one partition
        output_data (list): The output data vectors for the corresponding partition in the input
        output_scaler (scikit-learn scaler): The output data scaler
    Returns:
        model (scikit-learn model): The trained model
        [mae, mse, r2] (list of floats): The mean absolute error, mean square error, and the R^2 of the model
    """
    model = KernelRidge(kernel=params["kernel"], degree=params["degree"], 
                        alpha=params["alpha"], gamma=params["gamma"])
    model.fit(input_data, output_data)
    # Computing the error of the model
    prediction = model.predict(input_data)
    prediction_scaled = output_scaler.inverse_transform(prediction)
    target_output_data = output_scaler.inverse_transform(output_data)
    # Mean absolute error
    mae = mean_absolute_error(target_output_data, prediction_scaled)
    # Mean square error
    mse = mean_squared_error(target_output_data, prediction_scaled)
    # R^2
    r2 = r2_score(target_output_data, prediction_scaled)
    print('Training R^2:', r2)
    print('Training MSE:', mse, '  MAE:', mae)
    print('===========================================')

    return model, [mae, mse, r2]

def train(params):
    """
    This function is the main function that uses the previously defined 
    function to generate and train the data
    Args:
        params (dictionary):
            prefix: The prefix used to save the files.
            path_to_input_mats: The full path to input matrices files.
            path_to_output_mats: The full path to output matrices files.
            path_to_trajectory_xyz_file: The full path to trajectory xyz file. This is required for computation of the overlap
                                         matrix and solving generalized Kohn-Sham equations. Also, it will be required for computation of the molecular orbitals 
                                         overlap and time-overlap matrices.
            path_to_sample_files: The full path to sample files 
            input_property: The input property which appears in the middle of the name of the input files. It can take values of `kohn_sham`, `density`, `overlap`, and `hamiltonian` for xTB calculations.
            output_property: The output property which appears in the middle of the name of the output files. It can take values of `kohn_sham`, `density`, `overlap`, and `hamiltonian` for xTB calculations. 
            kernel: The KRR method kernels: `linear`, `poly` for polynomial kernels, and `rbf` for radial basis function kernel.
            degree: The degree of the kernel function in case of a `poly` kernel.
            alpha: This parameter represents the regularization strength in ridge regression. It controls the complexity of the model: a higher alpha value increases 
                   the amount of regularization, which in turn reduces the model's variance but might increase its bias. Specifically, it multiplies the identity matrix 
                   that is added to the kernel matrix in the ridge regression formula. This parameter helps to prevent overfitting by ensuring that the coefficients do not grow too large.
            gamma: This parameter is specific to certain types of kernels such as the radial basis function, sigmoid, or polynomial kernels. 
                   Gamma defines how much influence a single training input has. The larger the gamma, the closer other inputs must be to be affected.
            scaler: The scaler to scale the input and data. It can take `standard_scaler` and `minmax_scaler` values. The `standard_scaler` is recommended.
            save_models: A boolean flag for saving the trained model 
            path_to_save_models: The full path to save the model
            save_ml_ham: A boolean flag for saving the predicted Hamiltonian matrices 
            save_ao_overlap: A boolean flag for saving the atomic orbital overlap matrices
            save_ml_mos: A boolean flag for saving the molecular orbitals eigenvalues and eigenvectors
            partitioning_method: The partitioning of the Hamiltonian matrix. It can take values of `equal` and `atomic`. In the 
                                 equal partitioning, the upper triangular part of the Hamiltonian matrix is partitioned into `npartition` equal segments. 
                                 A separate model maps each input partition to its corresponding partition in the output matrix.
                                 The atomic partitioning method partitions the matrix based on the atomic angular momentum components of the matrix for each basis set.
            npartition: The number of partition in the `equal` partitioning method.
            memory_efficient: A boolean flag for memory efficiency of the calculations. This will remove the raw data and will remove them 
                              from the memory after they are processed.
            nprocs: The number of processors for computing the overlap matrices.
            write_wfn_file: A boolean flag for writing `wfn` files readable by CP2K.
            path_to_save_wfn_files: The full path to save the `wfn` files.
            is_periodic: The flag for periodicity of he system for computing the overlaps.
            A_cell_vector: A list containing the A cell vector.
            B_cell_vector: A list containing the B cell vector.
            C_cell_vector: A list containing the C cell vector.
            periodicity_type: The periodicity type for each direction of the periodic cell.
            lowest_orbital: The lowest orbital index to be saved. STARTS FROM 1
            highest_orbital: The highest orbital index to be saved. STARTS FROM 1
            res_dir: The full path to save the overlap, time-overlap, and energies 
            #TODO:
            train_parallel: A boolean flag for training models in parallel

    Returns: 
        models (list of models): A list of all trained models for each partition
        models_error (list of list of floats): Error data for each model including MAE, MSE, and R^2
        input_scalers (list of scalers): The list of all input scalers
        output_scalers (list of scalers): The list of all output scalers
    """
    critical_params = ['path_to_input_mats', 'path_to_output_mats', 'path_to_trajectory_xyz_file', 'path_to_sample_files']
    default_params = {'prefix': 'libra', 'input_property': 'kohn_sham', 'output_property': 'kohn_sham',
                      'kernel': 'linear', 'degree': 1, 'alpha': 1.0, 'gamma': 1.0, 'scaler': 'standard_scaler',
                      'save_models': True, 'path_to_save_models': './models', 'save_ml_ham': False, 'save_ao_overlap': False,
                      'save_ml_mos': False, 'partitioning_method': 'equal', 'npartition': 30, 'memory_efficient': True,
                      'nprocs': 2, 'write_wfn_file': False, 'path_to_save_wfn_files': './wfn_files', 'is_periodic': False,
                      'A_cell_vector': [25.0,0.0,0.0], 'B_cell_vector': [0.0,25.0,0.0], 'C_cell_vector': [0.0,0.0,25.0],
                      'periodicity_type': 'XYZ', 'lowest_orbital': 1, 'highest_orbital': 5, 'res_dir': './res', 'train_parallel': False
                     }
    # First load the default parameters etc
    comn.check_input(params, default_params, critical_params)
   
    # Find the training indices
    params["train_indices"] = find_indices(params)
    # Read the outputs 
    raw_output = read_data(params, params["path_to_output_mats"], params["prefix"]+"_ref_"+params["output_property"])
    # Read the inputs
    raw_input = read_data(params, params["path_to_input_mats"], params["prefix"]+"_"+params["input_property"])
    # Partition inputs
    params["input_partition"] = True # This auxiliary param is for correct partitioning using atomic partitioning
    partitioned_input = partition_data(params, raw_input)
    #print(len(partitioned_input), len(partitioned_input[0]), len(partitioned_input[0][0]))
    # Partition outputs
    params["input_partition"] = False
    partitioned_output = partition_data(params, raw_output)
    # Scale input data
    input_scalers, input_scaled = scale_data(params, partitioned_input)
    #print(len(input_scalers))
    # Scale output data
    output_scalers, output_scaled = scale_data(params, partitioned_output)
    # Train for each partition
    if params["train_parallel"]:
        raise("TODO!")
    # All models will be in this list
    else:
        models = []
        models_error = []
        for i in range(len(input_scaled)):
            model, model_error = train_partition(params, input_scaled[i], output_scaled[i], output_scalers[i])
            models.append(model)
            models_error.append(model_error)
        models_error = np.array(models_error)
    print("Average MAE for all models:", np.average(models_error[:,0]))
    print("Average MSE for all models:", np.average(models_error[:,1]))
    print("Average R^2 for all models:", np.average(models_error[:,2]))
    print("Done with training!!")
    if params["save_models"]:
        os.system(f"mkdir {params['path_to_save_models']}")
        np.save(f"{params['path_to_save_models']}/{params['prefix']}_models_error.npy", models_error)
        for i in range(len(models)):
            joblib.dump(models[i], F"{params['path_to_save_models']}/{params['prefix']}_model_{i}.joblib") 
        for i in range(len(output_scalers)):
            joblib.dump(output_scalers[i], F"{params['path_to_save_models']}/{params['prefix']}_output_scaler_{i}.joblib")
            joblib.dump(input_scalers[i], F"{params['path_to_save_models']}/{params['prefix']}_input_scaler_{i}.joblib")
    # Remove the input and output data to reduce the memory
    # del raw_input, raw_output, partitioned_input, partitioned_output
    # We also need the scalers when we want to use them for new data
    return models, models_error, input_scalers, output_scalers

def load_models(path_to_models, params):
    """
    This function loads the models that are already trained
    and return them with the input and output scalers as well for each partition.
    Args:
        path_to_models (string): The full path to saved models
        params (dict): The details of the parameters in the dictionary are given in the train function.
    Returns:
        models (list of models): A list of all trained models for each partition
        models_error (list of list of floats): Error data for each model including MAE, MSE, and R^2
        input_scalers (list of scalers): The list of all input scalers
        output_scalers (list of scalers): The list of all output scalers
    """
    models_error = np.load(f"{path_to_models}/{params['prefix']}_models_error.npy")
    models = []
    input_scalers = []
    output_scalers = []
    print("Loading models...")
    try:
        for i in range(params["npartition"]):
            models.append(joblib.load(f"{path_to_models}/{params['prefix']}_model_{i}.joblib"))
            input_scalers.append(joblib.load(f"{path_to_models}/{params['prefix']}_input_scaler_{i}.joblib"))
            output_scalers.append(joblib.load(f"{path_to_models}/{params['prefix']}_output_scaler_{i}.joblib"))
    except:
        raise("Could not load models!...")
    print("Done with loading models!...")
    return models, models_error, input_scalers, output_scalers



def compute_atomic_orbital_overlap_matrix(params, step):
    """
    This function computes the atomic orbital overlap
    matrix which will then be used to solve the generalized
    Kohn-Sham equations and compute the eigenvalues and
    eigenvectors. This is doen in the basis of the output data.
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        step (int): The step in the molecular dynamics trjectory
    Returns:
        AO_S (nparray): The atomic orbital overlap matrix
    """
    params["sample_molden_file"] = glob.glob(f"{params['path_to_sample_files']}/*ref*molden")[0]
    sample_molden_file = params["sample_molden_file"]
    path_to_trajectory = params["path_to_trajectory_xyz_file"]
    nprocs = params["nprocs"]
    is_periodic = params["is_periodic"]
    is_spherical = True
    molden_file_1 = F'temp_{step}.molden'
    molden_methods.write_molden_file(molden_file_1, sample_molden_file, path_to_trajectory, step)
    shell_1, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_1, is_spherical)
    AO_S = compute_overlaps(shell_1,shell_1,nprocs)
    if is_periodic:
        cell = []
        cell.append(params['A_cell_vector'])
        cell.append(params['B_cell_vector'])
        cell.append(params['C_cell_vector'])
        cell = np.array(cell)*units.Angst
        translational_vectors = params["translational_vectors"]
        for i1 in range(len(translational_vectors)):
            translational_vector = np.array(translational_vectors[i1])
            shell_1p, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_1, is_spherical, is_periodic, cell, translational_vector)
            AO_S += compute_overlaps(shell_1,shell_1p, nprocs)
    AO_S = data_conv.MATRIX2nparray(AO_S)
    os.system(F'rm  {molden_file_1}')
    new_indices = CP2K_methods.resort_ao_matrices(l_vals)
    #new_indices = molden_methods.resort_eigenvectors(l_vals)
    return AO_S[:,new_indices][new_indices,:]

def compute_mo_overlaps(params, eigenvectors_1, eigenvectors_2, step_1, step_2):
    """
    This function computes the overlap between two eigenvectors 
    where their geometries are given by step_1 and step_2 in 
    a molecular dynamics trajectory. This function can be used either 
    for computing the molecular orbital overlap matrix for one geometry
    or the time-overlap matrix of the molecular orbitals of two different
    geometries.
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        eigenvectors_1 (nparray): The first set of eigenvectors
        eigenvectors_2 (nparray): The second set of eigenvectors
        step_1 (int): The corresponding step in the molecular dynamics trajectory of eigenvectors_1
        step_2 (int): The corresponding step in the molecular dynamics trajectory of eigenvectors_2
    Returns:
        MO_overlap (nparray): The molecular orbital overlap matrix
    """
    params["sample_molden_file"] = glob.glob(f"{params['path_to_sample_files']}/*ref*molden")[0]
    sample_molden_file = params["sample_molden_file"]
    path_to_trajectory = params["path_to_trajectory_xyz_file"]
    nprocs = params["nprocs"]
    is_periodic = params["is_periodic"]
    is_spherical = True
    molden_file_1 = F'temp_{step_1}.molden'
    molden_methods.write_molden_file(molden_file_1, sample_molden_file, path_to_trajectory, step_1)
    shell_1, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_1, is_spherical)
    molden_file_2 = F'temp_{step_2}.molden'
    molden_methods.write_molden_file(molden_file_2, sample_molden_file, path_to_trajectory, step_2)
    shell_2, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_2, is_spherical)
    AO_S = compute_overlaps(shell_1,shell_2,nprocs)
    if is_periodic:
        cell = []
        cell.append(params['A_cell_vector'])
        cell.append(params['B_cell_vector'])
        cell.append(params['C_cell_vector'])
        cell = np.array(cell)*units.Angst
        translational_vectors = params["translational_vectors"]
        for i1 in range(len(translational_vectors)):
            translational_vector = np.array(translational_vectors[i1])
            shell_2p, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_2, is_spherical, is_periodic, 
                                                                          cell, translational_vector)
            AO_S += compute_overlaps(shell_1,shell_2p, nprocs)
    AO_S = data_conv.MATRIX2nparray(AO_S)
    os.system(F'rm {molden_file_1}')
    if molden_file_1!=molden_file_2:
        os.system(F'rm {molden_file_2}')
    new_indices = CP2K_methods.resort_ao_matrices(l_vals)
    AO_S = AO_S[:,new_indices][new_indices,:]
    MO_overlap = np.linalg.multi_dot([eigenvectors_1, AO_S, eigenvectors_2.T])
    #print(np.diag(MO_overlap))
    return MO_overlap

def find_indices_inputs(params):
    """
    The same as function find_indices_outputs, this function
    also finds the indices of input matrices
    This will gives us the initial step which will be a useful quantity 
    in indexing e.g. indexing the correct geometry in the molecular
    dynamics trajectory
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
    Returns:
        indices (list): The completed steps in the input directory
    """
    # We first find all the input files in the path_to_input_mats directory
    files = glob.glob(f'{params["path_to_input_mats"]}/*{params["input_property"]}*npy')
    indices = []
    for file in files:
        # Generate the index in a file
        index = file.split('/')[-1].split('_')[-1].replace('.npy','')
        indices.append(int(index))

    # Now sort the indices since they will be used
    indices = np.sort(indices)
    return list(indices)


def rebuild_matrix_from_partitions(params, partitions, output_shape):
    """
    This function is one of the most important here. It will
    rebuild a matix from its partitions by figuring out how it
    was originally partitioned.
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        partitions (list): The list of partitions of the upper triangular part of the Hamiltonian matrix
        output_shape (tuple): The shape of the output KS Ham matrix
    Returns:
        matrix (nparray): The rebuilt KS Ham matrix form partitions
    """
    # =========== 1- 'equal' or 'atomic' partitioning
    if params["partitioning_method"]=="equal" or params["partitioning_method"]=="atomic":
        upper_vector = np.concatenate(partitions, axis=0)
        #print(upper_vector.shape)
        # generate the upper indices of the matrix
        upper_indices = np.triu_indices(output_shape[0])
        matrix = np.zeros(output_shape)
        matrix[upper_indices] = upper_vector
        matrix = matrix + matrix.T - np.diag(matrix.diagonal())
    # =========== 2- 'block' partitioning
    elif params["partitioning_method"]=="block":
        print("To be implemented")

        
    return matrix


def compute_properties(params, models, input_scalers, output_scalers):
    """
    This function computes the molecular orbitals and their energy levels
    from the machine learned mapped Hamiltonian for each geometry.
    It also computes the overlap and time-overlap matrices between
    the computed molecular orbitals.
    It can output CP2K readable binary wfn files which contain the 
    basis set data and the molecular orbital coefficients, energies etc.
    This function does not return anything but it writes output results 
    to the params["res_dir"] which can then be used for computing the 
    excited states basis.
    Args:
        params (dict): The details of the parameters in the dictionary are given in the train function.
        models (list): The list that contains the trained models
        input_scalers (list): The list of all input scalers
        output_scalers (list): The list of all output scalers
    Returns:
        None
    """
    # First find the indices of the input data
    indices = find_indices_inputs(params)
    params["istep"] = indices[0]
    lowest_orbital = params["lowest_orbital"]
    highest_orbital = params["highest_orbital"]
    if params["write_wfn_file"]:
        try:
            params["sample_wfn_file"] = glob.glob(f"{params['path_to_sample_files']}/*ref*wfn")[0]
            basis_data, spin_data, eigen_vals_and_occ_nums, mo_coeffs = CP2K_methods.read_wfn_file(params["sample_wfn_file"])
            os.system(f'mkdir {params["path_to_save_wfn_files"]}')
        except:
            print("Sample wfn file not found or could not be read!")
            print("Will continue without writing the wfn files")
            params["write_wfn_file"] = False
    # Now we loop over the indices which are already sorted
    os.system(f'mkdir {params["res_dir"]}')
    for i, step in enumerate(indices):
        print("======================== \n Performing calculations for step ", step)
        input_mat = np.load(f'{params["path_to_input_mats"]}/{params["prefix"]}_{params["input_property"]}_{step}.npy')
        if i==0: 
            output_mat = np.load(f'{params["path_to_output_mats"]}/{params["prefix"]}_ref_{params["output_property"]}_{step}.npy')
        partitioned_input = partition_matrix(params, input_mat)
        #print(np.array(partitioned_input[0]).shape)
        #raise('  ')
        # Now apply the models to each partition
        outputs = []
        for j in range(len(input_scalers)):
            #print(input_scaled[j].reshape(1,-1).shape)
            input_scaled = input_scalers[j].transform(np.array(partitioned_input[j]).reshape(1,-1))
            #tmp = input_scaled.reshape(1,-1)
            #print(tmp.shape)
            output_scaled = models[j].predict(input_scaled)#.reshape(1,-1))
            output = output_scalers[j].inverse_transform(output_scaled)
            outputs.append(np.squeeze(output))
        #print(np.array(outputs).shape)
        ks_ham_mat = rebuild_matrix_from_partitions(params, outputs, output_mat.shape)
        #print(np.diag(tmp))
        #print('------')
        #print(np.diag(ks_ham_mat))
        #raise('   ')
        atomic_overlap = compute_atomic_orbital_overlap_matrix(params, step)
        eigenvalues, eigenvectors = CP2K_methods.compute_energies_coeffs(ks_ham_mat, atomic_overlap)
        if params["save_ml_ham"]:
            if not os.path.exists("ml_hams"):
                os.system(f"mkdir ml_hams")
            np.save(f"ml_hams/Ham_{step}.npy", ks_ham_mat)
        if params["save_ao_overlap"]:
            if not os.path.exists("ao_overlaps"):
                os.system(f"mkdir ao_overlaps")
            np.save(f"ao_overlaps/AO_S_{step}.npy", atomic_overlap)
        if params["save_ml_mos"]:
            if not os.path.exists("ml_mos"):
                os.system(f"mkdir ml_mos")
            np.save(f"ml_mos/eigenvectors_{step}.npy", eigenvectors) 
            np.save(f"ml_mos/eigenvalues_{step}.npy", eigenvalues) 
        # Compute the overlaps
        mo_overlap = compute_mo_overlaps(params, eigenvectors, eigenvectors, step, 
                                         step)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        # Compute the time-overlaps if we are at the next step
        zero_mat = np.zeros(mo_overlap.shape)
        if i>0:
            mo_time_overlap = compute_mo_overlaps(params, eigenvectors_prev, eigenvectors, step, step+1)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
            mo_time_overlap = data_conv.form_block_matrix(mo_time_overlap, zero_mat, zero_mat, mo_time_overlap)
            mo_time_overlap_sparse = sp.csc_matrix(mo_time_overlap)
            sp.save_npz(params["res_dir"]+F'/St_ks_{step-1}.npz', mo_time_overlap_sparse)
        energy_mat = np.diag(eigenvalues)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        #print(np.diag(energy_mat)*27.211385)
        
        energy_mat = data_conv.form_block_matrix(energy_mat, zero_mat, zero_mat, energy_mat)
        #output_name = F'ml_c60_b3lyp_size_{size}_{step+istep}-RESTART.wfn' 
        if params["write_wfn_file"]:
            output_name = F'{params["path_to_save_wfn_files"]}/ml_{params["prefix"]}-{step}-RESTART.wfn'
            CP2K_methods.write_wfn_file(output_name, basis_data, spin_data, eigen_vals_and_occ_nums, [eigenvectors[:,:]])
        eigenvectors_prev = eigenvectors
        # Writing the energy, overlap, and time-overlap files as in .npz file for other steps
        energy_mat_sparse = sp.csc_matrix(energy_mat)
        sp.save_npz(params["res_dir"]+F'/E_ks_{step}.npz', energy_mat_sparse)
        mo_overlap_sparse = sp.csc_matrix(mo_overlap)
        sp.save_npz(params["res_dir"]+F'/S_ks_{step}.npz', mo_overlap_sparse)


