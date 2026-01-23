# *********************************************************************************
# * Copyright (C) 2026 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: ci
   :platform: Unix, Windows
   :synopsis: this module implements functions for computing ci-level time-overlaps
.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

import numpy as np
from . import interfaces
#from .import 

def overlap(st_mo, data1, data2, params):
    """
    Compute the CI-state overlap matrix between two electronic-structure
    datasets using molecular-orbital time overlaps.

    This routine:
      1. Builds a common Slater-determinant basis from excited-state
         configurations of both datasets
      2. Constructs CI coefficient matrices in that common basis
      3. Computes SD and CSF overlap matrices (singlet)
      4. Projects the overlaps into the CI-state representation

    Parameters
    ----------
    st_mo : sparse matrix or array-like, shape (2*norb, 2*norb)
        Molecular-orbital time-overlap matrix in the spin–orbital basis.

    data1, data2 : tuple or list
        Electronic-structure data containers with the following layout::

            dataX[1] : list of list
                State-resolved configuration lists for excited states.
                dataX[1][i] contains configurations for excited state i+1.

            dataX[2] : list of list
                Corresponding CI amplitudes.
                dataX[2][i][j] is the amplitude of configuration j in
                excited state i+1.

        The ground state is assumed to be a pure reference determinant
        and is not included explicitly.

    params : dict
        Dictionary of required parameters:
            homo_indx : int
                HOMO index (1-based spatial orbital index)
            nocc : int
                Number of occupied orbitals below HOMO included in the window
            nvirt : int
                Number of virtual orbitals above HOMO included in the window
            nelec : int
                Total number of electrons
            nstates : int
                Total number of electronic states, including the ground state

    Returns
    -------
    numpy.ndarray
        CI-state overlap matrix of shape `(nstates, nstates)`.

    Notes
    -----
    - This routine is restricted to closed-shell singlet states.
    - Only singly excited configurations are assumed.
    - The CI overlap is computed as::

          S_CI = C1ᵀ · S_CSF · C2

      where `C1` and `C2` are CI coefficient matrices in a common CSF basis.
    """
    # ------------------------------------------------------------------
    # Extract and validate parameters
    # ------------------------------------------------------------------
    required_keys = {"homo_indx", "nocc", "nvirt", "nelec", "nstates"}
    missing = required_keys - params.keys()
    if missing:
        raise KeyError(f"Missing required parameters: {missing}")

    homo_indx = params["homo_indx"]
    nocc = params["nocc"]
    nvirt = params["nvirt"]
    nelec = params["nelec"]
    nstates = params["nstates"]

    if nstates <= 1:
        raise ValueError("nstates must include at least one excited state")

    # Orbital window (1-based spatial indices)
    lowest_orbital = homo_indx - nocc
    highest_orbital = homo_indx + nvirt

    # ------------------------------------------------------------------
    # Build common SD basis from excited states only
    # ------------------------------------------------------------------
    n_excited_states = nstates - 1

    common_sd_basis = interfaces.unique_confs(
        data1[1], data2[1], n_excited_states
    )

    # ------------------------------------------------------------------
    # CI coefficient matrices in the common SD basis
    # ------------------------------------------------------------------
    C1 = interfaces.ci_amplitudes_mtx(
        nstates, common_sd_basis, data1[1], data1[2]
    )
    C2 = interfaces.ci_amplitudes_mtx(
        nstates, common_sd_basis, data2[1], data2[2]
    )

    # ------------------------------------------------------------------
    # SD and CSF overlap matrices (singlet)
    # ------------------------------------------------------------------
    csf_ovlp, sd_ovlp = interfaces.sd_and_csf_overlaps_singlet(
        st_mo,
        lowest_orbital,
        highest_orbital,
        nelec,
        homo_indx,
        common_sd_basis,
    )

    # ------------------------------------------------------------------
    # CI-state overlap matrix
    # ------------------------------------------------------------------
    st_ci = C1.T @ csf_ovlp @ C2

    return st_ci

