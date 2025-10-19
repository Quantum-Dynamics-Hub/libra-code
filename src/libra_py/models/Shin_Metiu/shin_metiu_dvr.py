import numpy as np
import h5py
from scipy.special import erf
from functools import reduce
from liblibra_core import *

# ===============================================================
# Core components
# ===============================================================

def kinetic_energy_matrix(N, dr):
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


def V_en_matrix_sm1(r_list, R, L_sm1, Rc_sm1):
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

def V_en_matrix_sm2(r_list, R, L_sm2, a_sm2, af_sm2):
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


def V_en_deri_matrix_sm1(r_list, R, L_sm1, Rc_sm1):
    r_list = np.asarray(r_list)
    dV = np.zeros_like(r_list)
    for sigma in [+1, -1]:
        R_shift = R + sigma * L_sm1 / 2.0
        dV += -1.0 / np.abs(R_shift)**2 * np.sign(R_shift)
    R_shift = R - r_list
    dV -= 2.0 / np.sqrt(np.pi) / Rc_sm1 * np.exp(-(R_shift/Rc_sm1)**2) / np.abs(R_shift) * np.sign(R_shift)
    dV += erf(np.abs(R_shift)/Rc_sm1) / np.abs(R_shift)**2 * np.sign(R_shift)
    return np.diag(dV)

def V_en_deri_matrix_sm2(r_list, R, L_sm2, a_sm2, af_sm2):
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
    mu = R - np.asarray(r_list)
    return np.diag(mu)


# ===============================================================
# DVR Calculation (returns obj)
# ===============================================================

def compute_model(q, params, full_id=None):
    """
    Compute single-point Shin-Metiu DVR data and return as an obj.
    Compatible with later dynamics/polariton codes.
    """

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
    T = kinetic_energy_matrix(N, dr)
    H = T + V

    eigvals, eigvecs = np.linalg.eigh(H)
    # Transformations to adiabatic representation
    # <i| dV/dR |j>
    dV_adia = reduce(np.dot, (eigvecs.conj().T, dV, eigvecs))
    # Nonadiabatic couplings
    # <i| d/dR |j>
    nac = np.zeros_like(dV_adia)
    for p in range(len(eigvals)):
        for q in range(len(eigvals)):
            if p != q:
                nac[p, q] = dV_adia[p, q] / (eigvals[q] - eigvals[p])
    # d<i|V|j>/dR
    H_deri = dV_adia + np.dot(np.diag(eigvals), nac) - np.dot(nac, np.diag(eigvals))

    ham_adia = CMATRIX(nstates, nstates)
    hvib_adi = CMATRIX(nstates, nstates)
    basis_transform = CMATRIX(nstates, nstates)
    time_overlap_adi = CMATRIX(nstates, nstates)
    d1ham_adia = CMATRIXList()
    d1ham_adia.append(CMATRIX(nstates, nstates))
    dc1_adia = CMATRIXList()
    dc1_adia.append(CMATRIX(nstates, nstates))

    for i in range(nstates):
        ham_adia.set(i, i, eigvals[i] + 0.0j)
        hvib_adi.set(i, i, eigvals[i] + 0.0j)
        d1ham_adia[0].set(i, i, H_deri[i,i] + 0.0j)
        basis_transform.set(i, i, 1.0 + 0.0j)
        time_overlap_adi.set(i, i, 1.0 + 0.0j)
        for j in range(nstates):
            dc1_adia[0].set(i, j, nac[i,j]+0.0j)
    
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