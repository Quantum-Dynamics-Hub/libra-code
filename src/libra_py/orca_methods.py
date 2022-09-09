import os
import numpy as np



def resort_molog_eigenvectors(l_vals, p_perm, d_perm):
    """
    This function returns the resotring indices for resoting the MOLog 
    eigenvectors according to this order:
    
    MOLog order (example for Cd atom):
    2s, 3s | 3py, 3pz, 3px | 4py, 4pz, 4px | 4d-2, 4d-1, 4d0, 4d+1, 4d+2 | 5d-2, 5d-1, 5d0, 5d+1, 5d+2 | ...
    However, the atomic orbital overla computed from the libint version of Psi4 is not ordered 
    as above. The ordering is like this:
    2s, 3s | 3pz, 3px, 3py | 4pz, 4px, 4py | 4d0, 4d+1, 4d-1, 4d+2, 4d-2 | 5d0, 5d+1, 5d-1, 5d+2, 5d-2 | ...
    Therefore, we need to resort the eigenvectors to be able to use the code properly.
    
    Args:
    
        l_vals (list): A list containing the angular momentum values for atoms
                       in the order of the MOLog files.
                       
    Returns:
    
        new_indices (numpy array): The new indices that needs to be used for reordering.
    
    """
    # new indices
    new_indices = []
    # setting up the counter
    c = 0
    # loop over all the angular momentum values
    print(p_perm, d_perm)
    for i in range(len(l_vals)):
        l_val = l_vals[i]
        # find the reorder indices needed for this l_val
        reordered_indices = index_reorder(l_val, p_perm, d_perm)
        for j in range(len(reordered_indices)):
            # now append it by plus the counter since
            # we aim to do it for all the eigenvectors
            # and l values
            new_indices.append(c+reordered_indices[j])
        # increase the counter with t
        c += len(reordered_indices)
    # Return the new indices
    return new_indices



def index_reorder(l_val, p_perm, d_perm):
    """
    This function returns the new ordering of the angular momentum value based on the 
    order used in Psi4. Detailed explanation was given for the resort_molog_eigenvectors function.
    
    Args:
    
        l_val (integer): The angular momentum value.
                       
    Returns:
    
        new_order (numpy array): The new order of the indices for the l_val.
    
    """

    # for s orbital
    if l_val == 0:
        new_order = [1]
    # for p orbital
    elif l_val == 1:
        new_order = p_perm #[2,3,1]
    # for d orbital
    elif l_val == 2:
        new_order = d_perm #[1,2,3,4,5]
    # for f orbital
    elif l_val == 3:
        new_order = [4,5,3,6,2,7,1]
    # for g orbital
    elif l_val == 4:
        new_order = [5,6,4,7,3,8,2,9,1]

    # The indeices
    return np.array(new_order)-1



