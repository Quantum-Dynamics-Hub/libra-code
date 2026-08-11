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
.. module:: pyscf.implementations.test_casscf
   :platform: Unix, Windows
   :synopsis: Smoke test for PySCF CASSCF backend via the sequencing interface.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ES_Request, ES_Result, ES_Strategy, MolecularGeometry

import numpy as np


# Geometry coordinates are given explicitly in Bohr.
geom1 = MolecularGeometry(
    atom_labels=['He', 'H'],
    coords_bohr=np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.46379],
    ], dtype=np.float64)
)

geom2 = MolecularGeometry(
    atom_labels=['He', 'H'],
    coords_bohr=np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.88973],
    ], dtype=np.float64)
)

# The backend and the geometry container both use Bohr internally.
casscf = CASSCF(
    norbcas=2,
    nelecas=2,
    nroots=3,
    basis="sto-3g",
    charge=1,
    unit="Bohr",
)

request = ES_Request(
    n_singlets=3,
    n_triplets=0,
    H_soc=False,
    gradient_state="all",
    hessian_state=None,
    time_overlap=True,
    nacv=False,
)

result1 = ES_Result()
result2 = ES_Result()

casscf.compute_result(geom1, request, result1)
print("Result 1 H_el:", result1.H_el)
print("Result 1 gradients:", result1.gradients)
print("Result 1 time_overlap:", result1.time_overlap)

prev_state = casscf.copy()
casscf.compute_result(geom2, request, result2, previous=prev_state)
print("Result 2 H_el:", result2.H_el)
print("Result 2 gradients:", result2.gradients)
print("Result 2 time_overlap:", result2.time_overlap)

expected_overlap = np.array([
    [0.99692767, -0.01881889, 0.00438941],
    [-0.00177315, 0.96330492, -0.02404977],
    [-0.00817692, 0.02435447, 0.92265228],
], dtype=np.float64)
np.testing.assert_allclose(result2.time_overlap, expected_overlap, atol=1e-6, rtol=1e-6)

assert isinstance(casscf, ES_Strategy)
