# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: pyscf.implementations.test_casscf_nacv
   :platform: Unix, Windows
   :synopsis: NACV smoke test for the PySCF CASSCF backend via the ES interface.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>
"""

from __future__ import annotations

import numpy as np

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ES_Request, ES_Result, MolecularGeometry


NSTATES = 2
DISTANCES_BOHR = [6.0, 10.0]

basis_dict = {"Li": "sto-3g", "F": "6-311+g*"}
cas_list = [4, 7, 11, 14, 17]

casscf = CASSCF(
    norbcas=5,
    nelecas=2,
    nroots=NSTATES,
    basis=basis_dict,
    unit="Bohr",
    charge=0,
    cas_list=cas_list,
)

request = ES_Request(
    n_singlets=NSTATES,
    n_triplets=0,
    H_soc=False,
    gradient_state=None,
    hessian_state=None,
    nacv=True,
    time_overlap=True,
)

previous_strategy = None
for step, distance in enumerate(DISTANCES_BOHR):
    print(f"\n========== distance {step}: {distance:.2f} Bohr ==========")
    geom = MolecularGeometry(
        atom_labels=("Li", "F"),
        coords_bohr=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]], dtype=np.float64),
    )

    result = ES_Result()
    casscf.compute_result(geom, request, result, previous=previous_strategy)

    assert result.H_el is not None
    assert result.H_el.shape == (NSTATES,)
    print("energies:", result.H_el)

    assert result.nac_vectors is not None
    assert result.nac_vectors.shape == (NSTATES, NSTATES, 2, 3)
    print("nac_vectors shape:", result.nac_vectors.shape)

    assert result.time_overlap is not None
    assert result.time_overlap.shape == (NSTATES, NSTATES)
    print("time_overlap:\n", result.time_overlap)

    previous_strategy = casscf.copy()
