import numpy as np
import h5py
from scipy.special import erf
from functools import reduce
from liblibra_core import *

"""
.. module:: models_shin_meitu and dvr calculations
   :platform: Unix, Windows
   :synopsis: This module implements the Shin-meitu model (SM1 and SM2)
.. moduleauthor:: Yuchen Wang
.. module co-author:: Mohammad Shakiba, Alexey V. Akimov, Norah M. Hoffmann

"""

# ===============================================================
# Core components
# This module defines the key matrix elements used in a grid-based 
# representation of a 1D quantum system, such as a model for 
# electron-nuclear interactions. The matrices include:
# - Kinetic energy operator
# - Electron-nucleus potential energy (smooth models)
# - Derivatives of the potential (for forces or gradients)
# - Dipole operator
# ===============================================================

# ---------------------------------------------------------------
# Kinetic Energy Matrix
# ---------------------------------------------------------------
def kinetic_energy_matrix(N, dr):
    """
    Construct the kinetic energy matrix in a sine DVR (Discrete Variable Representation) basis.

    Parameters:
    - N : int - Number of grid points
    - dr : float - Grid spacing

    Returns:
    - T : ndarray (N x N) - Kinetic energy matrix
    """
    T = np.zeros((N, N))
    prefactor = 0.5
    for i in range(N):
        for j in range(N):
            if i == j:
                T[i, j] = prefactor * (np.pi ** 2) / (3 * dr ** 2) * (1 + 2 / N ** 2)
            else:
                delta = j - i
                sin_term = np.sin(np.pi * delta / N)
                T[i, j] = prefactor * 2 * (-1) ** delta * (np.pi ** 2) / ((dr * N * sin_term) ** 2)
    return T

# ---------------------------------------------------------------
# Model 1: Electron-Nucleus Potential Energy Matrix
# ---------------------------------------------------------------
def V_en_matrix_sm1(r_list, R, L_sm1, Rc_sm1):
    """
    Calculate the smoothened electron-nucleus potential energy matrix (model 1).

    Parameters:
    - r_list : array - Grid points
    - R : float - Nuclear position
    - L_sm1 : float - Softening length
    - Rc_sm1 : float - Smoothing parameter for erf

    Returns:
    - V : ndarray (N x N) - Diagonal potential matrix
    """
    r_list = np.asarray(r_list)
    V = np.zeros_like(r_list)
    for sigma in [+1, -1]:
        R_shift = R + sigma * L_sm1 / 2.0
        r_shift = r_list + sigma * L_sm1 / 2.0
        term1 = 1.0 / np.abs(R_shift)
        term2 = erf(np.abs(r_shift) / Rc_sm1) / np.abs(r_shift)
        V += term1 - term2
    V -= erf(np.abs(R - r_list) / Rc_sm1) / np.abs(R - r_list)
    return np.diag(V)

# ---------------------------------------------------------------
# Model 2: Electron-Nucleus Potential Energy Matrix
# ---------------------------------------------------------------
def V_en_matrix_sm2(r_list, R, L_sm2, a_sm2, af_sm2):
    """
    Calculate the smoothened electron-nucleus potential energy matrix (model 2).
    This version uses position-dependent smoothing parameters.

    Parameters:
    - r_list : array - Grid points
    - R : float - Nuclear position
    - L_sm2 : float - Softening length
    - a_sm2 : dict - {'+1': value, '-1': value} for different shifts
    - af_sm2 : float - Smoothing parameter for central term

    Returns:
    - V_matrix : ndarray (N x N) - Diagonal potential matrix
    """
    r_list = np.asarray(r_list)
    V = np.zeros_like(r_list)

    for sigma in [+1, -1]:
        R_shift = R + sigma * L_sm2 / 2.0
        r_shift = r_list + sigma * L_sm2 / 2.0
        term1 = 1.0 / np.abs(R_shift)
        term2 = erf(np.abs(r_shift) / a_sm2[sigma]) / np.abs(r_shift)
        V += term1 - term2

    V -= erf(np.abs(R - r_list) / af_sm2) / np.abs(R - r_list)

    V_matrix = np.diag(V)
    return V_matrix

# ---------------------------------------------------------------
# Derivative of Potential: Model 1
# ---------------------------------------------------------------
def V_en_deri_matrix_sm1(r_list, R, L_sm1, Rc_sm1):
    """
    Compute the derivative of the potential energy (model 1) with respect to R.

    Parameters:
    - r_list : array - Grid points
    - R : float - Nuclear position
    - L_sm1 : float - Softening length
    - Rc_sm1 : float - Smoothing parameter

    Returns:
    - dV : ndarray (N x N) - Diagonal matrix of potential derivatives
    """
    r_list = np.asarray(r_list)
    dV = np.zeros_like(r_list)
    for sigma in [+1, -1]:
        R_shift = R + sigma * L_sm1 / 2.0
        dV += -1.0 / np.abs(R_shift)**2 * np.sign(R_shift)
    R_shift = R - r_list
    dV -= 2.0 / np.sqrt(np.pi) / Rc_sm1 * np.exp(-(R_shift/Rc_sm1)**2) / np.abs(R_shift) * np.sign(R_shift)
    dV += erf(np.abs(R_shift)/Rc_sm1) / np.abs(R_shift)**2 * np.sign(R_shift)
    return np.diag(dV)

# ---------------------------------------------------------------
# Derivative of Potential: Model 2
# ---------------------------------------------------------------
def V_en_deri_matrix_sm2(r_list, R, L_sm2, a_sm2, af_sm2):
    """
    Compute the derivative of the potential energy (model 2) with respect to R.

    Parameters:
    - r_list : array - Grid points
    - R : float - Nuclear position
    - L_sm2 : float - Softening length
    - a_sm2 : dict - {'+1': value, '-1': value}
    - af_sm2 : float - Smoothing parameter for central term

    Returns:
    - dV_matrix : ndarray (N x N) - Diagonal matrix of potential derivatives
    """
    r_list = np.asarray(r_list)
    dV = np.zeros_like(r_list)

    for sigma in [+1, -1]:
        R_shift = R + sigma * L_sm2 / 2.0
        dV += -1.0 / np.abs(R_shift)**2 * np.sign(R_shift)

    R_shift = R - r_list
    dV -= 2.0 / np.sqrt(np.pi) / af_sm2 * np.exp(-(R_shift/af_sm2)**2) / np.abs(R_shift) * np.sign(R_shift)
    dV += erf(abs(R_shift)/af_sm2) / np.abs(R_shift)**2 * np.sign(R_shift)
    
    dV_matrix = np.diag(dV)
    return dV_matrix

def dipole_matrix(r_list, R):
    """
    Compute the dipole matrix (μ = R - r), diagonal in the position basis.

    Parameters:
    - r_list : array - Grid points
    - R : float - Nuclear position

    Returns:
    - mu : ndarray (N x N) - Diagonal dipole matrix
    """
    mu = R - np.asarray(r_list)
    return np.diag(mu)


# ===============================================================
# DVR Calculation (returns obj)
# This function computes the adiabatic Hamiltonian, its derivative with 
# respect to nuclear coordinate R, and nonadiabatic couplings.
# 
# The result is packaged into an object compatible with later 
# dynamics codes (e.g., surface hopping or polariton dynamics).
# ===============================================================

def compute_model(q, params, full_id=None):
    """
    Compute single-point Shin-Metiu DVR data and return as an object.

    Parameters:
    - q : float or object with get() method
        The current nuclear coordinate R (or an object providing it).
    - params : dict
        Dictionary of model parameters, including:
            - N : int, number of grid points
            - r_min, r_max : float, grid range
            - model : int (1 or 2), chooses between potential models
            - nstates : int, number of adiabatic states to keep
            - model-specific parameters (L_sm1, Rc_sm1, etc.)
    - full_id : optional
        Identifier used to retrieve correct nuclear coordinate from q.

    Returns:
    - obj : object
        Contains DVR results including adiabatic Hamiltonian,
        nonadiabatic couplings, and transformation matrices.
    """
    # --- Extract grid parameters --- 
    N = params["N"]
    nstates = params["nstates"]
    model = params["model"]
    r_grid = np.linspace(params["r_min"], params["r_max"], params["N"])
    dr = (params["r_max"] - params["r_min"]) / (params["N"] - 1)

    indx = 0
    if full_id is not None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]
    if isinstance(q, float):
        R = q
    else:
        R = q.get(0, indx)
    
    # --- Get nuclear coordinate R ---
    # Construct Hamiltonian
    if model == 1:
        L_sm1 = params["L_sm1"]
        Rc_sm1 = params["Rc_sm1"]
        V = V_en_matrix_sm1(r_grid, R, L_sm1, Rc_sm1)
        dV = V_en_deri_matrix_sm1(r_grid, R, L_sm1, Rc_sm1)
    elif model == 2:
        a_sm2 = params["a_sm2"]
        L_sm2 = params["L_sm2"]
        af_sm2 = params["af_sm2"]
        V = V_en_matrix_sm2(r_grid, R, L_sm2, a_sm2, af_sm2)
        dV = V_en_deri_matrix_sm2(r_grid, R, L_sm2, a_sm2, af_sm2)
    else:
        raise ValueError(F"Unknown model parameter {model} for Shin-Metiu Model.")
    # --- Construct Hamiltonian H = T + V --- 
    T = kinetic_energy_matrix(N, dr)
    H = T + V

    # --- Diagonalize H to get adiabatic states and energies ---
    eigvals, eigvecs = np.linalg.eigh(H)

    # --- Transformations to adiabatic representation --- 
    dV_adia = reduce(np.dot, (eigvecs.conj().T, dV, eigvecs))

    # --- Nonadiabatic couplings  <i| d/dR |j>  ---
    nac = np.zeros_like(dV_adia)
    for p in range(len(eigvals)):
        for q in range(len(eigvals)):
            if p != q:
                nac[p, q] = dV_adia[p, q] / (eigvals[q] - eigvals[p])
               
    # --- Compute dH/dR in adiabatic representation ---
    H_deri = dV_adia + np.dot(np.diag(eigvals), nac) - np.dot(nac, np.diag(eigvals))

    # --- Initialize output objects (custom complex matrix containers) ---
    ham_adia = CMATRIX(nstates, nstates)              # Adiabatic Hamiltonian
    hvib_adi = CMATRIX(nstates, nstates)              # Vibronic Hamiltonian (same as above)
    basis_transform = CMATRIX(nstates, nstates)       # Transformation matrix (identity here)
    time_overlap_adi = CMATRIX(nstates, nstates)      # Time overlap matrix (identity here)
    d1ham_adia = CMATRIXList()                        # First derivative of adiabatic Hamiltonian
    dc1_adi = CMATRIXList()                           # Nonadiabatic coupling vectors
    d1ham_adia.append(CMATRIX(nstates, nstates))
    dc1_adi.append(CMATRIX(nstates, nstates))

    for i in range(nstates):
        ham_adia.set(i, i, eigvals[i] + 0.0j)
        hvib_adi.set(i, i, eigvals[i] + 0.0j)
        d1ham_adia[0].set(i, i, H_deri[i,i] + 0.0j)
        basis_transform.set(i, i, 1.0 + 0.0j)
        time_overlap_adi.set(i, i, 1.0 + 0.0j)
        for j in range(nstates):
            dc1_adia[0].set(i, j, nac[i,j]+0.0j)
    
    # --- Create and return result object ---
    class tmp:
        pass
    obj = tmp()
    # Fill the obj
    obj.ham_adi = ham_adia
    obj.d1ham_adi = d1ham_adia
    obj.dc1_adi = dc1_adia
    obj.hvib_adi = hvib_adi
    obj.basis_transform = basis_transform
    obj.time_overlap_adi = time_overlap_adi

    return obj
