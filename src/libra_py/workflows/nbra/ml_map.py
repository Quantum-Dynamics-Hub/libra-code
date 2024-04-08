
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



def find_indices(params):
    """
    This function finds the indices of a set of file in a directory
    sorts them and return them in a list
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
    """
    data = []
    for i in params["train_indices"]:
        tmp = np.load(F'{path}/{prefix}_{params["input_property"]}_{i}.npy')
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
    3- 'atomwise': The section of a matrix which shows the interaction of each 
                   atom with all other atoms
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
    # =========== 3- 'atomwise' partitioning
    elif params["partitioning_method"]=="atomwise":
        print("To be implemented")

    return partitioned_matrix

def partition_data(params, data):
    """
    This function uses the partition_matrix to partition a set of data
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
    return (scalers, np.array(scaled_data))

def train_partition(params, input_data, output_data, output_scaler):
    """
    This function builds and trains a KRR model
    using input_data and output_data
    We will use the scaled input and output data to train the model
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
    """
    # Find the training indices
    params["train_indices"] = find_indices(params)
    # Read the outputs 
    raw_output = read_data(params, params["path_to_output_mats"], params["prefix"]+"_ref")
    # Read the inputs
    raw_input = read_data(params, params["path_to_input_mats"], params["prefix"])
    # Partition inputs
    partitioned_input = partition_data(params, raw_input)
    #print(len(partitioned_input), len(partitioned_input[0]), len(partitioned_input[0][0]))
    # Partition outputs
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

def train_parallel(params):
    """
    This function is the main function that uses the previously defined 
    function to generate and train the data but uses a parallel approach for 
    training each model for each partition
    TODO: This function speed is slower than the serial one! :)
    Need to create an auxiliary function for that and then use pool.map
    """
    # Find the training indices
    params["train_indices"] = find_indices(params)
    # Read the outputs 
    raw_output = read_data(params, params["path_to_output_mats"], params["prefix"])
    # Read the inputs
    raw_input = read_data(params, params["path_to_input_mats"], params["prefix"]+"_ref")
    # Partition inputs
    partitioned_input = partition_data(raw_input)
    #print(len(partitioned_input), len(partitioned_input[0]), len(partitioned_input[0][0]))
    # Partition outputs
    partitioned_output = partition_data(raw_output)
    # Scale input data
    input_scalers, input_scaled = scale_data(partitioned_input)
    #print(len(input_scalers))
    # Scale output data
    output_scalers, output_scaled = scale_data(partitioned_output)
    if params["memory_efficient"]:
        del raw_input, raw_output, partitioned_input, partitioned_output
    # Train for each partition
    # All models will be in this list
    """
    The only different part is here where we first 
    create a pool of processors and then use 
    pool.starmap to do the calculations
    We first create the set of arguments
    """
    arguments = []
    for i in range(params["npartition"]):
        arguments.append( (params, input_scaled[i], output_scaled[i], output_scalers[i]) )
    
    #with mp.Pool(processes=params["nprocs"]) as pool:
    #    results = pool.starmap(train_partition, arguments)
    pool = mp.Pool(processes=params["nprocs"])
    print('Started pool!')
    pool.starmap( train_partition, arguments )
    pool.close()
    pool.join()

    # To be completed
    print(results[0])


def compute_atomic_orbital_overlap_matrix(params, step):
    """
    This function computes the atomic orbital overlap
    matrix which will then be used to solve the generalized
    Kohn-Sham equations and compute the eigenvalues and
    eigenvectors. This is doen in the basis of the output data.
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
    """
    # =========== 1- 'equal' partitioning
    if params["partitioning_method"]=="equal":
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
    # =========== 3- 'atomwise' partitioning
    elif params["partitioning_method"]=="atomwise":
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
        for j in range(params["npartition"]):
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


