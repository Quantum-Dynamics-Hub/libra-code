# *********************************************************************************
# * Copyright (C) 2025 Yuchen Wang, Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: models_shin_meitu model with polaritonic effects
   :platform: Unix, Windows
   :synopsis: This module implements the polaritonic effects using Shin-metiu potential
   :Calculate the polaritonic states with the 2-state and 4-state models
    in reference paper.
    See (34)-(40) in J. Chem. Phys. 157, 104118 (2022)
.. moduleauthor:: Yuchen Wang
.. module co-author:: Mohammad Shakiba, Alexey V. Akimov, Norah M. Hoffmann

"""

import numpy as np
import h5py
from functools import reduce
from liblibra_core import *
import util.libutil as comn

# constants
from dataclasses import dataclass

HARTREE2EV = 27.2116  # constant stays as a constant


@dataclass
class DVRData:
    """
    Description of data: 
    * d_X = [\nabla X] and X_deri = \nabla [X], see reference paper

      solve cavity-free electronic DVR
    * eigvals: adiabatic electronic energies
    * eigvecs: adiabatic electronic wavefunctions
    * d_V: matrix of potential energy gradient [\nabla V]
    *  nac: non-adiabatic coupling matrix
    * mu: dipole moment matrix
    * d_mu: dipole moment gradient matrix [\nabla \mu]
    * mu_deri: derivative of dipole moment matrix elements \nabla [\mu]
    """
    R_grid: np.ndarray
    eigvals: np.ndarray
    d_V: np.ndarray
    nac: np.ndarray
    mu: np.ndarray
    d_mu: np.ndarray
    mu_deri: np.ndarray


def load_dvr_from_file(filename: str) -> DVRData:
    """
    Loads the DVR Shin-Meitu surfacefrom fitted potential file with HDF5 format

    Args:
        * filename (str): the name of the ".h5" file that contains data

    Returns: 
        DVRData : the digital form of the data stored in the file
    """

    with h5py.File(filename, 'r') as f:
        return DVRData(
            R_grid=f['R_grid'][:],
            eigvals=f['eigvals'][:],
            d_V=f['d_V'][:],
            nac=f['nac'][:],
            mu=f['mu'][:],
            d_mu=f['d_mu'][:],
            mu_deri=f['mu_deri'][:]
        )

def interpolateDVR_from_file(R: float, dvr_data: DVRData):
    """
    Interpolates grid-based DVR quantities to a specified nuclear coordinate.

    This function takes discrete DVR (Discrete Variable Representation) data
    loaded from an HDF5 file and performs 1D linear interpolation to evaluate
    electronic structure properties at an arbitrary position `R` along the
    nuclear coordinate grid.

    Args:
        R (float):
            Nuclear coordinate at which the properties should be evaluated.
        dvr_data (DVRData):
            An object containing DVR-computed quantities on a predefined
            coordinate grid. Must have attributes:
                - R_grid (1D array): Grid of nuclear coordinates.
                - eigvals (2D array [n_R, n_states]):
                    Electronic state energies.
                - d_V (3D array [n_R, n_states, n_states]):
                    Derivatives of the potential energy matrix.
                - nac (3D array [n_R, n_states, n_states]):
                    Nonadiabatic couplings.
                - mu (3D array [n_R, n_states, n_states]):
                    Dipole moment matrix elements.
                - d_mu (3D array [n_R, n_states, n_states]):
                    Derivatives of dipole moments.
                - mu_deri (3D array [n_R, n_states, n_states]):
                    Higher-order derivatives of dipole moments.

    Returns:
        tuple:
            A tuple containing:
            - eigvals (1D ndarray [n_states]):
                Interpolated electronic state energies.
            - d_V (2D ndarray [n_states, n_states]):
                Interpolated derivatives of the potential energy matrix.
            - nac (2D ndarray [n_states, n_states]):
                Interpolated nonadiabatic couplings.
            - mu (2D ndarray [n_states, n_states]):
                Interpolated dipole moment matrix.
            - d_mu (2D ndarray [n_states, n_states]):
                Interpolated dipole moment derivatives.
            - mu_deri (2D ndarray [n_states, n_states]):
                Interpolated higher-order dipole moment derivatives.

    Notes:
        - Interpolation is done independently for each matrix element using
          `numpy.interp` (1D linear interpolation).
        - This function assumes that the `R_grid` values are sorted in ascending
          order.
    """
    nelec_dim_dump = dvr_data.eigvals.shape[-1]

    eigvals = np.zeros((nelec_dim_dump))
    d_V = np.zeros((nelec_dim_dump, nelec_dim_dump))
    nac = np.zeros((nelec_dim_dump, nelec_dim_dump))
    mu = np.zeros((nelec_dim_dump, nelec_dim_dump))
    d_mu = np.zeros((nelec_dim_dump, nelec_dim_dump))
    mu_deri = np.zeros((nelec_dim_dump, nelec_dim_dump))

    for i in range(nelec_dim_dump):
        eigvals[i] = np.interp(R, dvr_data.R_grid, dvr_data.eigvals[:, i])
        for j in range(nelec_dim_dump):
            d_V[i, j] = np.interp(R, dvr_data.R_grid, dvr_data.d_V[:, i, j])
            nac[i, j] = np.interp(R, dvr_data.R_grid, dvr_data.nac[:, i, j])
            mu[i, j] = np.interp(R, dvr_data.R_grid, dvr_data.mu[:, i, j])
            d_mu[i, j] = np.interp(R, dvr_data.R_grid, dvr_data.d_mu[:, i, j])
            mu_deri[i, j] = np.interp(R, dvr_data.R_grid, dvr_data.mu_deri[:, i, j])

    return eigvals, d_V, nac, mu, d_mu, mu_deri

def polariton_info(R, dvr_data: DVRData, model='4-state',
                   g_c=0.005, omega_c=0.1, epsilon=1.0,
                   force_subspace=True, ndim_elec=None, ndim_ph=None):
    """
    Constructs polaritonic Hamiltonian and its derivatives in the adiabatic basis.

    This function interpolates electronic structure quantities from DVR data at a given
    nuclear coordinate `R`, then builds the coupled light–matter Hamiltonian for different
    model spaces ('2-state', '4-state', or 'general'). It also computes the corresponding
    derivatives with respect to `R` and the nonadiabatic coupling matrices in the adiabatic basis.

    Args:
        R (float):
            Nuclear coordinate at which the Hamiltonian should be evaluated.
        dvr_data (DVRData):
            DVR-computed quantities, with attributes:
                - R_grid (1D array)
                - eigvals, d_V, nac, mu, d_mu, mu_deri (grid-based arrays)
        model (str, optional):
            Model space choice: '2-state', '4-state', or 'general'. Default is '4-state'.
        g_c (float, optional):
            Light–matter coupling strength.
        omega_c (float, optional):
            Cavity photon frequency.
        epsilon (float, optional):
            Polarization factor (projection of the field on the dipole).
        force_subspace (bool, optional):
            Reserved for future use (restrict dynamics to a model subspace).
        ndim_elec (int, optional):
            Number of electronic states in the model. For '2-state' and '4-state' this is set automatically.
        ndim_ph (int, optional):
            Number of photon Fock states for 'general' model. Required if model == 'general'.

    Returns:
        tuple:
            - V_adia_fock (ndarray): Polaritonic Hamiltonian in adiabatic–Fock basis.
            - d_V_adia_fock (ndarray): Derivative of the polaritonic Hamiltonian with respect to `R`.
            - nac_adia_fock (ndarray): Nonadiabatic coupling matrix in adiabatic–Fock basis.

    Raises:
        AssertionError:
            If `model` is not one of the supported strings.
        NotImplementedError:
            If `model == 'general'` (not yet implemented).
        ValueError:
            If 'general' model is chosen without specifying `ndim_ph`.

    Notes:
        - The equations used here correspond to (27), (34), (36), (37), (38), and (39)
          in the reference paper.
        - For '2-state' and '4-state' models, photon zero-point energy shifts of 0.5*omega_c
          are included to match the reference.
    """

    assert model in ['2-state', '4-state', 'general'], "Model must be either '2-state' or '4-state'."

    if model == 'general':
        assert ndim_ph is not None and ndim_ph > 0

    if ndim_elec is None:
        ndim_elec = 2 if model in ['2-state', '4-state'] else ndim_elec
    # Load interpolated values
    eigvals, d_V, nac, mu, d_mu, mu_deri = interpolateDVR_from_file(R, dvr_data)

    # H_deri = d_H + np.dot(H, nac) - np.dot(nac, H)
    H_deri = d_V + np.dot(np.diag(eigvals), nac) - np.dot(nac, np.diag(eigvals))  # (27)

    # electronic part of D^2
    # assuming polarization direction aligned with dipole
    # (31d)
    D_square = epsilon**2 *g_c**2 / omega_c * np.dot(mu, mu)
    # gradient of D^2
    D_square_deri = epsilon**2 *g_c**2 / omega_c * (np.dot(mu, mu_deri) + np.dot(mu_deri, mu))

    # d_D_square = epsilon**2 *g_c**2 / omega_c * (np.dot(mu, d_mu) + np.dot(d_mu, mu))
    # D_square_deri = d_D_square + np.dot(D_square, nac) - np.dot(nac, D_square) # (27)

    # construct Hamiltonian in adiabatic-Fock basis
    # although the reference denote it as V but it is actually the Hamiltonian
    # Attention: for 2- and 4-state models, the Hamiltonian is shifted by 0.5*omega_c
    # to match the reference paper.
    V_adia_fock, d_V_adia_fock, nac_adia_fock = None, None, None
    if model == '2-state':
        # |e_0>, |g_1>
        # Eqn. (34)
        V_adia_fock = np.zeros((2, 2))
        V_adia_fock[0,0] = eigvals[1] + D_square[1,1]
        V_adia_fock[0,1] = g_c * epsilon * mu[1,0]
        V_adia_fock[1,0] = g_c * epsilon * mu[0,1]
        V_adia_fock[1,1] = eigvals[0] + D_square[0,0] + omega_c
        # Eqn. (36), in this case, d_V = H_deri
        V_deri_adia_fock = np.zeros((2, 2))
        V_deri_adia_fock[0,0] = H_deri[1,1] + D_square_deri[1,1]
        V_deri_adia_fock[0,1] = g_c * epsilon * mu_deri[1,0]
        V_deri_adia_fock[1,0] = g_c * epsilon * mu_deri[0,1]
        V_deri_adia_fock[1,1] = H_deri[0,0] + D_square_deri[0,0]
        d_V_adia_fock = V_deri_adia_fock.copy()

        nac_adia_fock = np.zeros((2, 2)) # TODO: need to define it
    elif model == '4-state':
        # |g_0>, |e_0>, |g_1>, |e_1>
        # Eqn. (37)
        V_adia_fock = np.zeros((4, 4))
        V_adia_fock[0,0] = eigvals[0] + D_square[0,0] + 0.5 * omega_c
        V_adia_fock[1,1] = eigvals[1] + D_square[1,1] + 0.5 * omega_c
        V_adia_fock[2,2] = eigvals[0] + D_square[0,0] + 1.5 * omega_c
        V_adia_fock[3,3] = eigvals[1] + D_square[1,1] + 1.5 * omega_c
        V_adia_fock[0,1] = D_square[0,1]
        V_adia_fock[2,3] = D_square[0,1]
        V_adia_fock[0,2] = g_c * epsilon * mu[0,0]
        V_adia_fock[0,3] = g_c * epsilon * mu[0,1]
        V_adia_fock[1,2] = g_c * epsilon * mu[1,0]
        V_adia_fock[1,3] = g_c * epsilon * mu[1,1]
        V_adia_fock[1,0] = V_adia_fock[0,1]
        V_adia_fock[3,2] = V_adia_fock[2,3]
        V_adia_fock[2,0] = V_adia_fock[0,2]
        V_adia_fock[3,0] = V_adia_fock[0,3]
        V_adia_fock[2,1] = V_adia_fock[1,2]
        V_adia_fock[3,1] = V_adia_fock[1,3]
        # Eqn. (38)
        nac_adia_fock = np.zeros((4, 4))
        nac_adia_fock[0,1] = nac[0,1]
        nac_adia_fock[2,3] = nac[0,1]
        nac_adia_fock[1,0] = nac[1,0]
        nac_adia_fock[3,2] = nac[1,0]
        assert np.linalg.norm(nac_adia_fock + nac_adia_fock.T) < 1e-10, \
            "Non-adiabatic coupling matrix is not anti-symmetric."
        # Eqn. (39)
        V_deri_adia_fock = np.zeros((4, 4))
        V_deri_adia_fock[0,0] = H_deri[0,0] + D_square_deri[0,0]
        V_deri_adia_fock[1,1] = H_deri[1,1] + D_square_deri[1,1]
        V_deri_adia_fock[2,2] = H_deri[0,0] + D_square_deri[0,0]
        V_deri_adia_fock[3,3] = H_deri[1,1] + D_square_deri[1,1]
        V_deri_adia_fock[0,1] = D_square_deri[0,1]
        V_deri_adia_fock[2,3] = D_square_deri[0,1]
        V_deri_adia_fock[0,2] = g_c * epsilon * mu_deri[0,0]
        V_deri_adia_fock[0,3] = g_c * epsilon * mu_deri[0,1]
        V_deri_adia_fock[1,2] = g_c * epsilon * mu_deri[1,0]
        V_deri_adia_fock[1,3] = g_c * epsilon * mu_deri[1,1]
        V_deri_adia_fock[1,0] = V_deri_adia_fock[0,1]
        V_deri_adia_fock[3,2] = V_deri_adia_fock[2,3]
        V_deri_adia_fock[2,0] = V_deri_adia_fock[0,2]
        V_deri_adia_fock[3,0] = V_deri_adia_fock[0,3]
        V_deri_adia_fock[2,1] = V_deri_adia_fock[1,2]
        V_deri_adia_fock[3,1] = V_deri_adia_fock[1,3]
        # (27)
        d_V_adia_fock = V_deri_adia_fock - np.dot(V_adia_fock, nac_adia_fock) + np.dot(nac_adia_fock, V_adia_fock)
    elif model == 'general':
        raise NotImplementedError("General model is not implemented yet.")
    else:
        raise ValueError("Model must be either '2-state', '4-state', or 'general'.")

    return V_adia_fock.copy(), d_V_adia_fock.copy(), nac_adia_fock.copy()  # eigvals_polariton, grad_polariton, nac_polariton, eigvecs_polariton


def compute_model(q, params, full_id):
    """
    Constructs a light–matter coupled Hamiltonian model for a given nuclear geometry.

    This function builds the adiabatic Hamiltonian, overlap matrix, and derivative
    couplings for a cavity–molecule system using precomputed DVR (Discrete Variable
    Representation) data and the polaritonic model specified by the user.

    Args:
        q:
            Coordinate accessor object with a `.get(coord_index, trajectory_index)`
            method for retrieving nuclear coordinates. `q.get(0, indx)` should
            return the current nuclear coordinate `R` for the system.
        params (dict):
            Dictionary of model parameters. Recognized keys:
                - "g_c" (float, optional):
                    Light–matter coupling strength. Default: 0.005
                - "omega_c" (float, optional):
                    Cavity photon frequency. Default: 0.1
                - "epsilon" (float, optional):
                    Polarization factor (projection of the field on the dipole).
                    Default: 1.0
                - "nstates" (int, optional):
                    Number of adiabatic states in the model (2 or 4). Default: 4
                - "model" (int):
                    Index selecting which DVR data file to use:
                        0 → "dvr_sm1.h5"
                        1 → "dvr_sm2.h5"
                    Any other value raises ValueError.
        full_id:
            Identifier for the current trajectory or computation. Passed to
            `Cpp2Py` to extract an index for coordinate retrieval.

    Returns:
        obj (object):
            A Python object with the following attributes:
                - ham_dia (CMATRIX):
                    Adiabatic Hamiltonian matrix (nstates × nstates).
                - ovlp_dia (CMATRIX):
                    Overlap matrix in adiabatic representation (identity matrix).
                - d1ham_dia (CMATRIXList):
                    List of first derivatives of the Hamiltonian with respect
                    to nuclear coordinates.
                - dc1_dia (CMATRIXList):
                    List of first-order derivative couplings
                    (empty here; initialized as zero matrices).

    Raises:
        ValueError:
            If `model` is not 0 or 1.
        KeyError:
            If required parameters are missing from `params`.

    Notes:
        - DVR data is loaded from an HDF5 file once per call (no global caching).
        - The actual polaritonic Hamiltonian is constructed by calling
          `polariton_info`, which returns the Hamiltonian, its derivative,
          and the nonadiabatic coupling matrix in the adiabatic basis.
        - Overlap matrix `S_adia` is set to the identity matrix in this implementation.
    """

    critical_params = []
    default_params = {"g_c": 0.005, "omega_c": 0.1, "epsilon": 1.0, "nstates": 4}
    comn.check_input(params, default_params, critical_params)

    g_c = params["g_c"]
    omega_c = params["omega_c"]
    epsilon = params["epsilon"]
    nstates = params["nstates"]
    model = params["model"]

    if model == 0:
        dvr_file = "dvr_sm1.h5"
    elif model == 1:
        dvr_file = "dvr_sm2.h5"
    else:
        raise ValueError(f"Unknown nstates {nstates}")

    model_polariton = '2-state' if nstates == 2 else '4-state'

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    R = q.get(0, indx)

    # Load DVR data once, no globals
    dvr_data = load_dvr_from_file(dvr_file)

    V_dia_fock, dV_dia_fock, nac_adia_fock = polariton_info(
        R, dvr_data, model_polariton, g_c, omega_c, epsilon
    )

    H_adia = CMATRIX(nstates, nstates)
    S_adia = CMATRIX(nstates, nstates)
    d1ham_adia = CMATRIXList()
    d1ham_adia.append(CMATRIX(nstates, nstates))
    dc1_adia = CMATRIXList()
    dc1_adia.append(CMATRIX(nstates, nstates))

    for i in range(nstates):
        S_adia.set(i, i, 1.0 + 0.0j)
        for j in range(nstates):
            H_adia.set(i, j, V_dia_fock[i,j]+0.0j)
            d1ham_adia[0].set(i, j, dV_dia_fock[i,j]+0.0j)


    class tmp:
        pass

    obj = tmp()
    obj.ham_dia = H_adia
    obj.ovlp_dia = S_adia
    obj.d1ham_dia = d1ham_adia
    obj.dc1_dia = dc1_adia 


    return obj

