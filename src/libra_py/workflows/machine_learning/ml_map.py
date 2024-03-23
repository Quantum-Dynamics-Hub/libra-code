"""
This module performs ML calculations by mapping to be completed....
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
import pickle 
import joblib
# import tensorflow as tf
# from tensorflow.keras import layers, models

from sklearn.preprocessing import StandardScaler
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_squared_error, accuracy_score, mean_absolute_error, r2_score
from liblibra_core import *
import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.packages.dftbplus.methods as DFTB_methods
from libra_py import units, molden_methods, data_conv



def load_data(filename):
    """
    """
    if ".npy" in filename:
        data = np.load(filename)
    elif ".npz" in filename:
        data = np.array(sp.load_npz(filename).todense().real)
    return data


def compute_atomic_orbital_overlap_matrix(params, step):
    """
    This 
    """
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
    This 
    """
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
            shell_2p, l_vals = molden_methods.molden_file_to_libint_shell(molden_file_2, is_spherical, is_periodic, cell, translational_vector)
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




# I guess this function's name is not general and because of this I keep it in this module 
# otherwise I would move it to libra_py.data_conv
def upper_vector_to_symmetric_nparray(upper_vector, upper_indices, mat_shape):
    """
    This function gets the upper triangular part of a matrix as a vector and retuns a symmetric matrix
    Args:
        upper_vector (nparray): The upper part of a matrix
        upper_indices (nparray): The indices of the upper part of the matrix
        mat_shape (tuple): The shape of the original numpy array
    Returns:
        matrix (nparray): The symmetric matix built based on the upper triangular matrix
    """
    matrix = np.zeros(mat_shape)
    matrix[upper_indices] = upper_vector
    matrix = matrix + matrix.T - np.diag(matrix.diagonal())
    return matrix


def run_ml_map(params):
    """
    This function runs calculations for ml
    """
    # ========================= Part 1: Extracting the data
    print("Loading input matrices...")
    input_mats = load_data(params["path_to_input_mat"])
    print("Loading output matrices...")
    output_mats = load_data(params["path_to_output_mat"])
    print("Done with loading inputs and outputs...")

    # ========================= Part 2: Preparation of the inputs and outputs for training
    # Creating random shuffling of indices
    shuffled_indices = np.arange(input_mats.shape[0])
    np.random.shuffle(shuffled_indices)
    print(shuffled_indices, shuffled_indices.shape)
    # Saving the shuffled indices
    np.save('shuffled_indices.npy', shuffled_indices)
    # Reading the block indices of the input and output matrices
    # In fact, here we read the data block related to two atoms
    input_sample_angular_momentum_file = params["input_sample_angular_momentum_file"]
    output_sample_angular_momentum_file = params["output_sample_angular_momentum_file"]
    try:
        input_block_indices = CP2K_methods.atom_components_cp2k(input_sample_angular_momentum_file)
        print("Angular momentum components for input matrices extracted from CP2K calculations...")
    except:
        input_block_indices = DFTB_methods.atom_components_dftb(input_sample_angular_momentum_file)
        print("Angular momentum components for input matrices extracted from DTFB+ calculations...")

    try:
        output_block_indices = CP2K_methods.atom_components_cp2k(output_sample_angular_momentum_file)
        print("Angular momentum components for output matrices extracted from CP2K calculations...")
    except:
        output_block_indices = DFTB_methods.atom_components_dftb(output_sample_angular_momentum_file)
        print("Angular momentum components for output matrices extracted from DTFB+ calculations...")

    # The size of the training set
    size = params["training_set_size"]
    train_indices = shuffled_indices[0:size]
    test_indices = shuffled_indices[size:]
    # Append the first and last geometries for better predictions (this is an interpolation model)
    if max(shuffled_indices) not in train_indices:
        np.append(train_indices, max(shuffled_indices))
    if min(shuffled_indices) not in train_indices:
        np.append(train_indices, min(shuffled_indices))

    # In this part we vectorize the inputs and outputs to feed to model
    # The off-diagonal blocks
    inputs_train = []
    inputs_test = []
    outputs_train = []
    outputs_test = []
    # The diagonal blocks
    inputs_train_diag = []
    inputs_test_diag = []
    outputs_train_diag = []
    outputs_test_diag = []
    # Off-diagonal blocks: inputs
    for i in range(len(input_block_indices)):
        for j in range(len(input_block_indices)):
            if j>i:
                s1 = input_block_indices[i][0]
                e1 = input_block_indices[i][-1]+1
                s2 = input_block_indices[j][0]
                e2 = input_block_indices[j][-1]+1
                
                tmp1 = input_mats[train_indices,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                inputs_train.append(tmp2)
                
                tmp1 = input_mats[test_indices,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                inputs_test.append(tmp2)
    # Off-diagonal blocks: outputs
    for i in range(len(output_block_indices)):
        for j in range(len(output_block_indices)):
            if j>i:
                s1 = output_block_indices[i][0]
                e1 = output_block_indices[i][-1]+1
                s2 = output_block_indices[j][0]
                e2 = output_block_indices[j][-1]+1
                
                tmp1 = output_mats[train_indices,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                outputs_train.append(tmp2)
                
                tmp1 = output_mats[test_indices,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                outputs_test.append(tmp2)
                
    print("Shape of the vectorized off-diagonal inputs and outputs:")
    print(len(inputs_train), len(inputs_train[0]), len(inputs_train[0][0]))
    print(len(outputs_train), len(outputs_train[0]), len(outputs_train[0][0]))
    print(len(inputs_test), len(inputs_test[0]), len(inputs_test[0][0]))
    print(len(outputs_test), len(outputs_test[0]), len(outputs_test[0][0]))

    # Diagonal blocks: inputs
    for i in range(len(input_block_indices)):
        s1 = input_block_indices[i][0]
        e1 = input_block_indices[i][-1]+1
        tmp = input_mats[0][input_block_indices[i],:][:,input_block_indices[i]]
        upper_indices = np.triu_indices(tmp.shape[0])
        
        tmp1 = input_mats[train_indices, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        inputs_train_diag.append(tmp2)
        
        tmp1 = input_mats[test_indices, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        inputs_test_diag.append(tmp2)
        
    # Diagonal blocks: outputs
    for i in range(len(output_block_indices)):
        s1 = output_block_indices[i][0]
        e1 = output_block_indices[i][-1]+1
        tmp = output_mats[0][input_block_indices[i],:][:,input_block_indices[i]]
        upper_indices = np.triu_indices(tmp.shape[0])
        
        tmp1 = output_mats[train_indices, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        outputs_train_diag.append(tmp2)
        
        tmp1 = output_mats[test_indices, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        outputs_test_diag.append(tmp2)
   
    print("Shape of the vectorized diagonal inputs and outputs:")
    print(len(inputs_train_diag), len(inputs_train_diag[0]), len(inputs_train_diag[0][0]))
    print(len(outputs_train_diag), len(outputs_train_diag[0]), len(outputs_train_diag[0][0]))
    print(len(inputs_test_diag), len(inputs_test_diag[0]), len(inputs_test_diag[0][0]))
    print(len(outputs_test_diag), len(outputs_test_diag[0]), len(outputs_test_diag[0][0]))
    
    # Now, we are ready to scale the values
    # The scalers for each block
    input_scalers = []
    output_scalers = []

    scaler = params["scaler"]

    inputs_train_scaled = []
    inputs_test_scaled = []
    outputs_train_scaled = []
    outputs_test_scaled = []
    for i in range(len(inputs_train)):
        if scaler.lower()=="standard_scaler":
            input_scaler = StandardScaler()
            output_scaler = StandardScaler()
        elif scaler.lower()=="min_max_scaler":
            input_scaler = MinMaxScaler()
            output_scaler = MinMaxScaler()
        else:
            raise("Scaler not recognized. Try standard_scaler or min_max_scaler.")
        input_scaler.fit(inputs_train[i])
        output_scaler.fit(outputs_train[i])
        
        # Scaling the training set
        inputs_train_scaled.append(input_scaler.transform(inputs_train[i]))
        outputs_train_scaled.append(output_scaler.transform(outputs_train[i]))
        # Scaling the testing set
        inputs_test_scaled.append(input_scaler.transform(inputs_test[i]))
        outputs_test_scaled.append(output_scaler.transform(outputs_test[i]))
        input_scalers.append(input_scaler)
        output_scalers.append(output_scaler)


    input_scalers_diag = []
    output_scalers_diag = []
    
    inputs_train_diag_scaled = []
    inputs_test_diag_scaled = []
    outputs_train_diag_scaled = []
    outputs_test_diag_scaled = []
    for i in range(len(inputs_train_diag)):
        if scaler.lower()=="standard_scaler":
            input_scaler = StandardScaler()
            output_scaler = StandardScaler()
        elif scaler.lower()=="min_max_scaler":
            input_scaler = MinMaxScaler()
            output_scaler = MinMaxScaler()
        else:
            raise("Scaler not recognized. Try standard_scaler or min_max_scaler.")
        input_scaler.fit(inputs_train_diag[i])
        output_scaler.fit(outputs_train_diag[i])
        
        # Scaling the training set
        inputs_train_diag_scaled.append(input_scaler.transform(inputs_train_diag[i]))
        outputs_train_diag_scaled.append(output_scaler.transform(outputs_train_diag[i]))
        # Scaling the testing set
        inputs_test_diag_scaled.append(input_scaler.transform(inputs_test_diag[i]))
        outputs_test_diag_scaled.append(output_scaler.transform(outputs_test_diag[i]))
        input_scalers_diag.append(input_scaler)
        output_scalers_diag.append(output_scaler)

    # ========================= Part 3: Training with training data
    # =========== Training off-diagonal blocks
    # All models for each block
    models = []
    # The error of each model for training and testing sets
    model_train_error = []
    model_test_error = []
    for i in range(len(inputs_train)):
        # Creating the model
        if params["method"].lower()=="krr":
            model = KernelRidge(kernel=params["kernel"], degree=params["degree"], alpha=params["alpha"], gamma=params["gamma"]) 
        else:
            raise("Only KRR method is implemented for now!")
        # Training the model
        model.fit(inputs_train_scaled[i], outputs_train_scaled[i])
        # Computing the error of the model
        # === Training set
        predictions_train = model.predict(inputs_train_scaled[i])
        predictions_train_scaled = output_scalers[i].inverse_transform(predictions_train)
        # Mean absolute error
        mae = mean_absolute_error(outputs_train[i], predictions_train_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_train[i], predictions_train_scaled)
        # R^2
        r2 = r2_score(outputs_train[i], predictions_train_scaled)
        model_train_error.append([mae, mse, r2])
        print(F'Off-diagonal block model {i} ---> Training R2:', r2)
        print(F'Off-diagonal block model {i} ---> Training    MSE:', mse, '  MAE:', mae)
        # === Testing set
        predictions_test = model.predict(inputs_test_scaled[i])
        predictions_test_scaled = output_scalers[i].inverse_transform(predictions_test)
        # Mean absolute error
        mae = mean_absolute_error(outputs_test[i], predictions_test_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_test[i], predictions_test_scaled)
        # R^2 
        r2 = r2_score(outputs_test[i], predictions_test_scaled)
        model_test_error.append([mae,mse,r2])
        print(F'Off-diagonal block model {i} ---> Testing R2:', r2)
        print(F'Off-diagonal block model {i} ---> Testing     MSE:', mse, '  MAE:', mae)
        print('====================================================')
        models.append(model)

    # =========== Training diagonal blocks
    # All models for each block
    models_diag = []
    # The error of each model for training and testing sets
    model_train_diag_error = []
    model_test_diag_error = []
    for i in range(len(inputs_train_diag)):
        # Creating the model
        if params["method"].lower()=="krr":
            model = KernelRidge(kernel=params["kernel"], degree=params["degree"], alpha=params["alpha"], gamma=params["gamma"]) 
        else:
            raise("Only KRR method is implemented for now!")
        # Training the model
        model.fit(inputs_train_diag_scaled[i], outputs_train_diag_scaled[i])
        # Computing the error of the model
        # === Training set
        predictions_train = model.predict(inputs_train_diag_scaled[i])
        predictions_train_scaled = output_scalers_diag[i].inverse_transform(predictions_train)
        # Mean absolute error
        mae = mean_absolute_error(outputs_train_diag[i], predictions_train_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_train_diag[i], predictions_train_scaled)
        # R^2
        r2 = r2_score(outputs_train_diag[i], predictions_train_scaled)
        model_train_diag_error.append([mae, mse, r2])
        print(F'Diagonal block model {i} ---> Training R2:', r2)
        print(F'Diagonal block model {i} ---> Training    MSE:', mse, '  MAE:', mae)
        # === Testing set
        predictions_test = model.predict(inputs_test_diag_scaled[i])
        predictions_test_scaled = output_scalers_diag[i].inverse_transform(predictions_test)
        # Mean absolute error
        mae = mean_absolute_error(outputs_test_diag[i], predictions_test_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_test_diag[i], predictions_test_scaled)
        # R^2 
        r2 = r2_score(outputs_test_diag[i], predictions_test_scaled)
        model_test_diag_error.append([mae,mse,r2])
        print(F'Diagonal block model {i} ---> Testing R2:', r2)
        print(F'Diagonal block model {i} ---> Testing     MSE:', mse, '  MAE:', mae)
        print('====================================================')
        models_diag.append(model)
    os.system(F"mkdir {params['path_to_save_model']}")
    os.system(F"mkdir {params['path_to_save_wfn_files']}")
    os.system(F"mkdir {params['res_dir']}")
    if params["save_model"]:
        # Saving the models
        for i in range(len(models)):
            joblib.dump(models[i], F'{params["path_to_save_model"]}/off_diag_model_{i}.pkl')
        for i in range(len(models_diag)):
            joblib.dump(models_diag[i], F'{params["path_to_save_model"]}/diag_model_{i}.pkl')


    # ========================= Part 4: Using the model
    if params["write_wfn_file"]:
        try:
            basis_data, spin_data, eigen_vals_and_occ_nums, mo_coeffs = CP2K_methods.read_wfn_file(params["sample_wfn_file"])
        except:
            raise("Sample wfn file not found or could not be read!")
    istep = params["istep"]
    lowest_orbital = params["lowest_orbital"]
    highest_orbital = params["highest_orbital"]
    project = params["project"]
    for step in range(0, len(input_mats)):
        print(F"Calculating properties for step {step+istep}")
        c = 0
        ks_mat = np.zeros(output_mats[step].shape)
        for i in range(len(input_block_indices)):
            for j in range(len(input_block_indices)):
                if j>i:
                    s1 = input_block_indices[i][0]
                    e1 = input_block_indices[i][-1]+1
                    s2 = input_block_indices[j][0]
                    e2 = input_block_indices[j][-1]+1
                    block = input_mats[step, s1:e1, s2:e2]
                    block_shape = block.shape
                    block_ = block.reshape(1,block_shape[0]*block_shape[1])
                    #print(block.shape)
                    block_scaled = input_scalers[c].transform(block_)
                    block_predict = models[c].predict(block_scaled)
                    block_predict = output_scalers[c].inverse_transform(block_predict)
                    block_predict = block_predict.reshape(block.shape)
                    s1 = output_block_indices[i][0]
                    e1 = output_block_indices[i][-1]+1
                    s2 = output_block_indices[j][0]
                    e2 = output_block_indices[j][-1]+1
                    #print(block_predict)
                    ks_mat[s1:e1, s2:e2] = block_predict
                    c += 1
        # print(ks_mat[0,:])
        ks_mat = ks_mat + ks_mat.T
        for i in range(len(input_block_indices)):
            s1 = input_block_indices[i][0]
            e1 = input_block_indices[i][-1]+1
            block = input_mats[step, s1:e1, s1:e1]
            upper_indices = np.triu_indices(block.shape[0])
            tmp1 = input_mats[train_indices, s1:e1, s1:e1]
            block_ = block[upper_indices[0], upper_indices[1]]
            #print(block.shape, block_.shape)
            #print(block.shape)
            block_shape = block_.shape
            block_ = block_.reshape(1,block_shape[0])
            block_scaled = input_scalers_diag[i].transform(block_)
            block_predict = models_diag[i].predict(block_scaled)
            block_predict = output_scalers_diag[i].inverse_transform(block_predict)
            block_predict = upper_vector_to_symmetric_nparray(block_predict, upper_indices, block.shape)
            s1 = output_block_indices[i][0]
            e1 = output_block_indices[i][-1]+1
            ks_mat[s1:e1, s1:e1] = block_predict
        # print(ks_mat)
        # print(converged_ks_mats[0])
        # ====== Computing the overlap matrix from the trajectory
        atomic_overlap = compute_atomic_orbital_overlap_matrix(params, step+istep)
        eigenvalues, eigenvectors = CP2K_methods.compute_energies_coeffs(ks_mat, atomic_overlap)
        mo_overlap = compute_mo_overlaps(params, eigenvectors, eigenvectors, step+istep, step+istep)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        # Making the double spin format
        zero_mat = np.zeros(mo_overlap.shape)
        mo_overlap = data_conv.form_block_matrix(mo_overlap, zero_mat, zero_mat, mo_overlap)
        if step>0:
            mo_time_overlap = compute_mo_overlaps(params, eigenvectors_prev, eigenvectors, step+istep, step+istep+1)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
            mo_time_overlap = data_conv.form_block_matrix(mo_time_overlap, zero_mat, zero_mat, mo_time_overlap)
            mo_time_overlap_sparse = sp.csc_matrix(mo_time_overlap)
            sp.save_npz(params["res_dir"]+F'/St_ks_{step-1+istep}.npz', mo_time_overlap_sparse)
        energy_mat = np.diag(eigenvalues)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        energy_mat = data_conv.form_block_matrix(energy_mat, zero_mat, zero_mat, energy_mat)
        #output_name = F'ml_c60_b3lyp_size_{size}_{step+istep}-RESTART.wfn' 
        if params["write_wfn_file"]:
            output_name = F'{params["path_to_save_wfn_files"]}/ml_{project}-{step+istep}-RESTART.wfn' 
            CP2K_methods.write_wfn_file(output_name, basis_data, spin_data, eigen_vals_and_occ_nums, [eigenvectors[:,:]])
        eigenvectors_prev = eigenvectors
        # Writing the energy, overlap, and time-overlap files as in .npz file for other steps
        energy_mat_sparse = sp.csc_matrix(energy_mat)
        sp.save_npz(params["res_dir"]+F'/E_ks_{step+istep}.npz', energy_mat_sparse)
        mo_overlap_sparse = sp.csc_matrix(mo_overlap)
        sp.save_npz(params["res_dir"]+F'/S_ks_{step+istep}.npz', mo_overlap_sparse)




def run_ml_map_v2(params):
    """
    This function runs calculations for ml
    """
    # ========================= Part 1: Extracting the data
    #print("Loading input matrices...")
    #input_mats = load_data(params["path_to_input_mat"])
    #print("Loading output matrices...")
    #output_mats = load_data(params["path_to_output_mat"])
    #print("Done with loading inputs and outputs...")
    train_indices = params["train_indices"]
    all_indices = params["all_indices"]
    print(train_indices)
    print(all_indices)
    train_inputs = []
    train_outputs = []
    all_inputs = []
    prefix = params["prefix"]
    input_property = params["input_property"]
    output_property = params["output_property"] 
    path_to_input_mats = params["path_to_input_mats"]
    path_to_output_mats = params["path_to_output_mats"]
    input_mats = []
    all_input_mats = []
    output_mats = []
    for i in train_indices:
        tmp = np.load(F'{path_to_input_mats}/{prefix}_{input_property}_{i}.npy')
        input_mats.append(tmp)
        tmp = np.load(F'{path_to_output_mats}/{prefix}_{output_property}_{i}.npy')
        output_mats.append(tmp)
    input_mats = np.array(input_mats)
    print(input_mats.shape)
    output_mats = np.array(output_mats)
    print(output_mats.shape)
    
    for i in all_indices:
        #if i not in train_indices:
        tmp = np.load(F'{path_to_input_mats}/{prefix}_{input_property}_{i}.npy')
        all_input_mats.append(tmp)
    all_input_mats = np.array(all_input_mats)
    print(all_input_mats.shape) 
    # ========================= Part 2: Preparation of the inputs and outputs for training
    # Creating random shuffling of indices
    #shuffled_indices = np.arange(input_mats.shape[0])
    #np.random.shuffle(shuffled_indices)
    #print(shuffled_indices, shuffled_indices.shape)
    # Saving the shuffled indices
    #np.save('shuffled_indices.npy', shuffled_indices)
    # Reading the block indices of the input and output matrices
    # In fact, here we read the data block related to two atoms
    input_sample_angular_momentum_file = params["input_sample_angular_momentum_file"]
    output_sample_angular_momentum_file = params["output_sample_angular_momentum_file"]
    try:
        input_block_indices = CP2K_methods.atom_components_cp2k(input_sample_angular_momentum_file)
        print("Angular momentum components for input matrices extracted from CP2K calculations...")
    except:
        input_block_indices = DFTB_methods.atom_components_dftb(input_sample_angular_momentum_file)
        print("Angular momentum components for input matrices extracted from DTFB+ calculations...")

    try:
        output_block_indices = CP2K_methods.atom_components_cp2k(output_sample_angular_momentum_file)
        print("Angular momentum components for output matrices extracted from CP2K calculations...")
    except:
        output_block_indices = DFTB_methods.atom_components_dftb(output_sample_angular_momentum_file)
        print("Angular momentum components for output matrices extracted from DTFB+ calculations...")

    # The size of the training set
    #size = params["training_set_size"]
    #train_indices = shuffled_indices[0:size]
    #test_indices = shuffled_indices[size:]
    # Append the first and last geometries for better predictions (this is an interpolation model)
    #if max(shuffled_indices) not in train_indices:
    #    np.append(train_indices, max(shuffled_indices))
    #if min(shuffled_indices) not in train_indices:
    #    np.append(train_indices, min(shuffled_indices))

    # In this part we vectorize the inputs and outputs to feed to model
    # The off-diagonal blocks
    inputs_train = []
    inputs_test = []
    outputs_train = []
    outputs_test = []
    # The diagonal blocks
    inputs_train_diag = []
    inputs_test_diag = []
    outputs_train_diag = []
    outputs_test_diag = []
    # Off-diagonal blocks: inputs
    for i in range(len(input_block_indices)):
        for j in range(len(input_block_indices)):
            if j>i:
                s1 = input_block_indices[i][0]
                e1 = input_block_indices[i][-1]+1
                s2 = input_block_indices[j][0]
                e2 = input_block_indices[j][-1]+1
                
                tmp1 = input_mats[:,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                inputs_train.append(tmp2)
                
                tmp1 = all_input_mats[:,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                inputs_test.append(tmp2)

    # Off-diagonal blocks: outputs
    for i in range(len(output_block_indices)):
        for j in range(len(output_block_indices)):
            if j>i:
                s1 = output_block_indices[i][0]
                e1 = output_block_indices[i][-1]+1
                s2 = output_block_indices[j][0]
                e2 = output_block_indices[j][-1]+1
                
                tmp1 = output_mats[:,s1:e1,s2:e2]
                tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                outputs_train.append(tmp2)
                
                #tmp1 = output_mats[test_indices,s1:e1,s2:e2]
                #tmp2 = np.reshape(tmp1, (tmp1.shape[0], -1))
                #outputs_test.append(tmp2)
                
    print("Shape of the vectorized off-diagonal inputs and outputs:")
    print(len(inputs_train), len(inputs_train[0]), len(inputs_train[0][0]))
    print(len(outputs_train), len(outputs_train[0]), len(outputs_train[0][0]))
    print(len(inputs_test), len(inputs_test[0]), len(inputs_test[0][0]))
    #print(len(outputs_test), len(outputs_test[0]), len(outputs_test[0][0]))

    # Diagonal blocks: inputs
    for i in range(len(input_block_indices)):
        s1 = input_block_indices[i][0]
        e1 = input_block_indices[i][-1]+1
        tmp = input_mats[0][input_block_indices[i],:][:,input_block_indices[i]]
        upper_indices = np.triu_indices(tmp.shape[0])
        
        tmp1 = input_mats[:, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        inputs_train_diag.append(tmp2)
        
        tmp1 = all_input_mats[:, s1:e1, s1:e1]
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        inputs_test_diag.append(tmp2)
        
    # Diagonal blocks: outputs
    for i in range(len(output_block_indices)):
        s1 = output_block_indices[i][0]
        e1 = output_block_indices[i][-1]+1
        tmp = output_mats[0][input_block_indices[i],:][:,input_block_indices[i]]
        upper_indices = np.triu_indices(tmp.shape[0])
        
        tmp1 = output_mats[:, s1:e1, s1:e1] 
        tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        outputs_train_diag.append(tmp2)
        
        #tmp1 = output_mats[test_indices, s1:e1, s1:e1]
        #tmp2 = tmp1[:, upper_indices[0], upper_indices[1]]
        #outputs_test_diag.append(tmp2)
   
    print("Shape of the vectorized diagonal inputs and outputs:")
    print(len(inputs_train_diag), len(inputs_train_diag[0]), len(inputs_train_diag[0][0]))
    print(len(outputs_train_diag), len(outputs_train_diag[0]), len(outputs_train_diag[0][0]))
    print(len(inputs_test_diag), len(inputs_test_diag[0]), len(inputs_test_diag[0][0]))
    #print(len(outputs_test_diag), len(outputs_test_diag[0]), len(outputs_test_diag[0][0]))
    
    # Now, we are ready to scale the values
    # The scalers for each block
    input_scalers = []
    output_scalers = []

    scaler = params["scaler"]

    inputs_train_scaled = []
    inputs_test_scaled = []
    outputs_train_scaled = []
    #outputs_test_scaled = []
    for i in range(len(inputs_train)):
        if scaler.lower()=="standard_scaler":
            input_scaler = StandardScaler()
            output_scaler = StandardScaler()
        elif scaler.lower()=="min_max_scaler":
            input_scaler = MinMaxScaler()
            output_scaler = MinMaxScaler()
        else:
            raise("Scaler not recognized. Try standard_scaler or min_max_scaler.")
        input_scaler.fit(inputs_train[i])
        output_scaler.fit(outputs_train[i])
        
        # Scaling the training set
        inputs_train_scaled.append(input_scaler.transform(inputs_train[i]))
        outputs_train_scaled.append(output_scaler.transform(outputs_train[i]))
        # Scaling the testing set
        inputs_test_scaled.append(input_scaler.transform(inputs_test[i]))
        #outputs_test_scaled.append(output_scaler.transform(outputs_test[i]))
        input_scalers.append(input_scaler)
        output_scalers.append(output_scaler)


    input_scalers_diag = []
    output_scalers_diag = []
    
    inputs_train_diag_scaled = []
    inputs_test_diag_scaled = []
    outputs_train_diag_scaled = []
    #outputs_test_diag_scaled = []
    for i in range(len(inputs_train_diag)):
        if scaler.lower()=="standard_scaler":
            input_scaler = StandardScaler()
            output_scaler = StandardScaler()
        elif scaler.lower()=="min_max_scaler":
            input_scaler = MinMaxScaler()
            output_scaler = MinMaxScaler()
        else:
            raise("Scaler not recognized. Try standard_scaler or min_max_scaler.")
        input_scaler.fit(inputs_train_diag[i])
        output_scaler.fit(outputs_train_diag[i])
        
        # Scaling the training set
        inputs_train_diag_scaled.append(input_scaler.transform(inputs_train_diag[i]))
        outputs_train_diag_scaled.append(output_scaler.transform(outputs_train_diag[i]))
        # Scaling the testing set
        inputs_test_diag_scaled.append(input_scaler.transform(inputs_test_diag[i]))
        #outputs_test_diag_scaled.append(output_scaler.transform(outputs_test_diag[i]))
        input_scalers_diag.append(input_scaler)
        output_scalers_diag.append(output_scaler)

    # ========================= Part 3: Training with training data
    # =========== Training off-diagonal blocks
    # All models for each block
    models = []
    # The error of each model for training and testing sets
    model_train_error = []
    #model_test_error = []
    for i in range(len(inputs_train)):
        # Creating the model
        if params["method"].lower()=="krr":
            model = KernelRidge(kernel=params["kernel"], degree=params["degree"], alpha=params["alpha"], gamma=params["gamma"]) 
        else:
            raise("Only KRR method is implemented for now!")
        # Training the model
        model.fit(inputs_train_scaled[i], outputs_train_scaled[i])
        # Computing the error of the model
        # === Training set
        predictions_train = model.predict(inputs_train_scaled[i])
        predictions_train_scaled = output_scalers[i].inverse_transform(predictions_train)
        # Mean absolute error
        mae = mean_absolute_error(outputs_train[i], predictions_train_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_train[i], predictions_train_scaled)
        # R^2
        r2 = r2_score(outputs_train[i], predictions_train_scaled)
        model_train_error.append([mae, mse, r2])
        print(F'Off-diagonal block model {i} ---> Training R2:', r2)
        print(F'Off-diagonal block model {i} ---> Training    MSE:', mse, '  MAE:', mae)
        # === Testing set
        #predictions_test = model.predict(inputs_test_scaled[i])
        #predictions_test_scaled = output_scalers[i].inverse_transform(predictions_test)
        ## Mean absolute error
        #mae = mean_absolute_error(outputs_test[i], predictions_test_scaled)
        ## Mean square error
        #mse = mean_squared_error(outputs_test[i], predictions_test_scaled)
        ## R^2 
        #r2 = r2_score(outputs_test[i], predictions_test_scaled)
        #model_test_error.append([mae,mse,r2])
        #print(F'Off-diagonal block model {i} ---> Testing R2:', r2)
        #print(F'Off-diagonal block model {i} ---> Testing     MSE:', mse, '  MAE:', mae)
        print('====================================================')
        models.append(model)

    # =========== Training diagonal blocks
    # All models for each block
    models_diag = []
    # The error of each model for training and testing sets
    model_train_diag_error = []
    #model_test_diag_error = []
    for i in range(len(inputs_train_diag)):
        # Creating the model
        if params["method"].lower()=="krr":
            model = KernelRidge(kernel=params["kernel"], degree=params["degree"], alpha=params["alpha"], gamma=params["gamma"]) 
        else:
            raise("Only KRR method is implemented for now!")
        # Training the model
        model.fit(inputs_train_diag_scaled[i], outputs_train_diag_scaled[i])
        # Computing the error of the model
        # === Training set
        predictions_train = model.predict(inputs_train_diag_scaled[i])
        predictions_train_scaled = output_scalers_diag[i].inverse_transform(predictions_train)
        # Mean absolute error
        mae = mean_absolute_error(outputs_train_diag[i], predictions_train_scaled)
        # Mean square error
        mse = mean_squared_error(outputs_train_diag[i], predictions_train_scaled)
        # R^2
        r2 = r2_score(outputs_train_diag[i], predictions_train_scaled)
        model_train_diag_error.append([mae, mse, r2])
        print(F'Diagonal block model {i} ---> Training R2:', r2)
        print(F'Diagonal block model {i} ---> Training    MSE:', mse, '  MAE:', mae)
        # === Testing set
        #predictions_test = model.predict(inputs_test_diag_scaled[i])
        #predictions_test_scaled = output_scalers_diag[i].inverse_transform(predictions_test)
        ## Mean absolute error
        #mae = mean_absolute_error(outputs_test_diag[i], predictions_test_scaled)
        ## Mean square error
        #mse = mean_squared_error(outputs_test_diag[i], predictions_test_scaled)
        ## R^2 
        #r2 = r2_score(outputs_test_diag[i], predictions_test_scaled)
        #model_test_diag_error.append([mae,mse,r2])
        #print(F'Diagonal block model {i} ---> Testing R2:', r2)
        #print(F'Diagonal block model {i} ---> Testing     MSE:', mse, '  MAE:', mae)
        #print('====================================================')
        models_diag.append(model)
    os.system(F"mkdir {params['path_to_save_model']}")
    os.system(F"mkdir {params['path_to_save_wfn_files']}")
    os.system(F"mkdir {params['res_dir']}")
    if params["save_model"]:
        # Saving the models
        for i in range(len(models)):
            joblib.dump(models[i], F'{params["path_to_save_model"]}/off_diag_model_{i}.pkl')
        for i in range(len(models_diag)):
            joblib.dump(models_diag[i], F'{params["path_to_save_model"]}/diag_model_{i}.pkl')


    # ========================= Part 4: Using the model
    if params["write_wfn_file"]:
        try:
            basis_data, spin_data, eigen_vals_and_occ_nums, mo_coeffs = CP2K_methods.read_wfn_file(params["sample_wfn_file"])
        except:
            raise("Sample wfn file not found or could not be read!")
    istep = params["istep"]
    lowest_orbital = params["lowest_orbital"]
    highest_orbital = params["highest_orbital"]
    project = params["project"]
    for step in range(0, len(all_input_mats)):
        print(F"Calculating properties for step {step+istep}")
        c = 0
        ks_mat = np.zeros(output_mats[0].shape)
        for i in range(len(input_block_indices)):
            for j in range(len(input_block_indices)):
                if j>i:
                    s1 = input_block_indices[i][0]
                    e1 = input_block_indices[i][-1]+1
                    s2 = input_block_indices[j][0]
                    e2 = input_block_indices[j][-1]+1
                    block = all_input_mats[step, s1:e1, s2:e2]
                    block_shape = block.shape
                    block_ = block.reshape(1,block_shape[0]*block_shape[1])
                    #print(block.shape)
                    block_scaled = input_scalers[c].transform(block_)
                    block_predict = models[c].predict(block_scaled)
                    block_predict = output_scalers[c].inverse_transform(block_predict)
                    block_predict = block_predict.reshape(block.shape)
                    s1 = output_block_indices[i][0]
                    e1 = output_block_indices[i][-1]+1
                    s2 = output_block_indices[j][0]
                    e2 = output_block_indices[j][-1]+1
                    #print(block_predict)
                    ks_mat[s1:e1, s2:e2] = block_predict
                    c += 1
        # print(ks_mat[0,:])
        ks_mat = ks_mat + ks_mat.T
        for i in range(len(input_block_indices)):
            s1 = input_block_indices[i][0]
            e1 = input_block_indices[i][-1]+1
            block = all_input_mats[step, s1:e1, s1:e1]
            upper_indices = np.triu_indices(block.shape[0])
            #tmp1 = all_input_mats[train_indices, s1:e1, s1:e1]
            block_ = block[upper_indices[0], upper_indices[1]]
            #print(block.shape, block_.shape)
            #print(block.shape)
            block_shape = block_.shape
            block_ = block_.reshape(1,block_shape[0])
            block_scaled = input_scalers_diag[i].transform(block_)
            block_predict = models_diag[i].predict(block_scaled)
            block_predict = output_scalers_diag[i].inverse_transform(block_predict)
            block_predict = upper_vector_to_symmetric_nparray(block_predict, upper_indices, block.shape)
            s1 = output_block_indices[i][0]
            e1 = output_block_indices[i][-1]+1
            ks_mat[s1:e1, s1:e1] = block_predict
        # print(ks_mat)
        # print(converged_ks_mats[0])
        # ====== Computing the overlap matrix from the trajectory
        atomic_overlap = compute_atomic_orbital_overlap_matrix(params, step+istep)
        eigenvalues, eigenvectors = CP2K_methods.compute_energies_coeffs(ks_mat, atomic_overlap)
        mo_overlap = compute_mo_overlaps(params, eigenvectors, eigenvectors, step+istep, step+istep)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        # Making the double spin format
        zero_mat = np.zeros(mo_overlap.shape)
        mo_overlap = data_conv.form_block_matrix(mo_overlap, zero_mat, zero_mat, mo_overlap)
        if step>0:
            mo_time_overlap = compute_mo_overlaps(params, eigenvectors_prev, eigenvectors, step+istep, step+istep+1)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
            mo_time_overlap = data_conv.form_block_matrix(mo_time_overlap, zero_mat, zero_mat, mo_time_overlap)
            mo_time_overlap_sparse = sp.csc_matrix(mo_time_overlap)
            sp.save_npz(params["res_dir"]+F'/St_ks_{step-1+istep}.npz', mo_time_overlap_sparse)
        energy_mat = np.diag(eigenvalues)[lowest_orbital-1:highest_orbital,:][:,lowest_orbital-1:highest_orbital]
        energy_mat = data_conv.form_block_matrix(energy_mat, zero_mat, zero_mat, energy_mat)
        #output_name = F'ml_c60_b3lyp_size_{size}_{step+istep}-RESTART.wfn' 
        if params["write_wfn_file"]:
            output_name = F'{params["path_to_save_wfn_files"]}/ml_{project}-{step+istep}-RESTART.wfn' 
            CP2K_methods.write_wfn_file(output_name, basis_data, spin_data, eigen_vals_and_occ_nums, [eigenvectors[:,:]])
        eigenvectors_prev = eigenvectors
        # Writing the energy, overlap, and time-overlap files as in .npz file for other steps
        energy_mat_sparse = sp.csc_matrix(energy_mat)
        sp.save_npz(params["res_dir"]+F'/E_ks_{step+istep}.npz', energy_mat_sparse)
        mo_overlap_sparse = sp.csc_matrix(mo_overlap)
        sp.save_npz(params["res_dir"]+F'/S_ks_{step+istep}.npz', mo_overlap_sparse)













