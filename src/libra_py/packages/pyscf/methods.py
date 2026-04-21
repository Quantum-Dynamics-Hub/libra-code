# *********************************************************************************
# * Copyright (C) 2026  Jieyang Gu and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: methods
   :platform: Unix, Windows
   :synopsis: This module implements the Libra/PySCF interface function

.. moduleauthor::
       Alexey V. Akimov, Jieyang Gu

"""


import os, sys, math, re, struct, copy, subprocess
import numpy as np
from liblibra_core import MATRIX, CMATRIX, CMATRIXList, Py2Cpp_int, Cpp2Py, Random

# Fisrt, we add the location of the library to test to the PYTHON path
from libra_py.packages.pyscf.implementations.cisd import CISD
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry

from libra_py.packages.cp2k import methods as cp2k
from libra_py import data_conv
from libra_py import units

class tmp:
    pass


def pyscf_compute_adi(q, params, full_id):

    # ================= Decode trajectory index =================
    Id = Cpp2Py(full_id)
    itraj = Id[-1]

    # ================= Extract coordinates =================
    coords = q.col(itraj)
    coordinates = data_conv.MATRIX2nparray(coords, float).reshape(-1,3)/units.Angst # in Angstrom

    ndof = coords.num_of_rows
    nat = ndof // 3

    # ================= Safe param access =================
    params.setdefault("is_first_time", {})
    params.setdefault("act_state", {})
    #params.setdefault("pyscf_obj", {})
    params.setdefault("coords_prev", {})

    is_first_time = params["is_first_time"].get(itraj, True)
    act_state = params["act_state"].get(itraj, 0)
    
    _basis = params.get("basis", "sto-3g")
    _charge = params.get("charge", 0)
    nstates = params.get("nstates", 2)
    method = params.get("method", "casscf")
    # The default is CAS(2,2) 
    _norbcas = params.get("norbcas", 2)
    _nelecas = params.get("nelecas", 2)

    # ================= Read parameters =================
    dt = float(params.get("dt", 41.0))
    _atom_labels = params["atom_labels"]
    
    # ================= Previous coordinates =================
    
    if is_first_time:
        coords_prev = copy.deepcopy(coordinates)
    else:
        coords_prev = copy.deepcopy(params["coords_prev"].get(itraj))

    geom_prev = MolecularGeometry(atom_labels = _atom_labels, coords_angstrom=np.array(coords_prev) )
    geom = MolecularGeometry(atom_labels = _atom_labels, coords_angstrom=np.array(coordinates) )
    
        
    pyscf_obj = None
    if method=="casscf":
        pyscf_obj = CASSCF(norbcas=_norbcas, nelecas=_nelecas, nroots=nstates, basis=_basis, charge=_charge)
    elif method=="cisd":
        pyscf_obj = CISD(nroots=nstates, basis=_basis, charge=_charge)     
    
    
    pyscf_obj.set_geom_and_run_hf(geom_prev)
    _ = [pyscf_obj.compute_energy(root) for root in range(nstates)]
    _ = [pyscf_obj.compute_gradient(root) for root in range(nstates)]

    pyscf_obj.set_geom_and_run_hf(geom)
    energies = [pyscf_obj.compute_energy(root) for root in range(nstates)]
    grad = [pyscf_obj.compute_gradient(root) for root in range(nstates)]

    # ================= Compute overlaps =================
    st_ci = pyscf_obj.time_overlap_matrix(nstates)

    print(F"coords_prev = {coords_prev}")
    print(F"coordinates = {coordinates}")

    # ================= Build object =================
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)

    # ================= Populate Hamiltonian =================
    for i in range(nstates):
        obj.ham_adi.set(i, i, energies[i] * (1.0+0.0j) )
        obj.hvib_adi.set(i, i, energies[i] * (1.0+0.0j) )
        obj.basis_transform.set(i, i, 1.0+0.0j)

        for j in range(nstates):
            obj.time_overlap_adi.set(i, j, float(st_ci[i, j]) * (1.0+0.0j) )

    # ================== Forces ===============================
    obj.d1ham_adi = CMATRIXList()
    for idof in range(ndof):
        obj.d1ham_adi.append(CMATRIX(nstates, nstates))

    for iatom in range(nat):
        for i in range(nstates):
            obj.d1ham_adi[3 * iatom + 0].set(i, i, grad[i][iatom, 0] * (1.0 + 0.0j))
            obj.d1ham_adi[3 * iatom + 1].set(i, i, grad[i][iatom, 1] * (1.0 + 0.0j))
            obj.d1ham_adi[3 * iatom + 2].set(i, i, grad[i][iatom, 2] * (1.0 + 0.0j))

    # ================= Compute derivative couplings =================
    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = ( obj.time_overlap_adi.get(i, j) - obj.time_overlap_adi.get(j, i) ) / (2.0 * dt)
            obj.hvib_adi.set(i, j, -1j * dij)
            obj.hvib_adi.set(j, i, +1j * dij)            
            
    # ================= Store state =================
    #params["pyscf_obj"][itraj] = copy.deepcopy(pyscf_obj)
    #params["pyscf_obj"][itraj] = pyscf_obj
    params["coords_prev"][itraj] = copy.deepcopy(coordinates)
    params["is_first_time"][itraj] = False

    return obj

