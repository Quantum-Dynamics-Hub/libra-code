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
   :synopsis: Smoke test for PySCF CASSCF backend.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

import sys
from pathlib import Path

# Allow running this script directly with proper package root in sys.path
#if __name__ == "__main__" and __package__ is None:
#    file_path = Path(__file__).resolve()
#    for parent in file_path.parents:
#        if parent.name == "src":
#            sys.path.insert(0, str(parent))
#            break
#    else:
#        raise RuntimeError("Could not locate src/ directory on path for libra_py import")

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry
import numpy as np

NTRAJ = 2
NSTATES = 3

geom_step0 = [
    MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7746]])),
    MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9000]])),
]

geom_step1 = [
    MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0000]])),
    MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.1000]])),
]

# instantiate the CASSCF strategy class with parameters for the test
casscf = CASSCF(norbcas=2, nelecas=2, nroots=NSTATES, basis='sto-3g', charge=1, ntraj=NTRAJ)

for traj_id, geom in enumerate(geom_step0):
    # First geometry for each trajectory initializes that trajectory's state slot.
    casscf.set_geom_and_run_hf(geom, traj_id=traj_id)

    energies = [casscf.compute_energy(root, traj_id=traj_id) for root in range(NSTATES)]
    print(f'Trajectory {traj_id} energies at step 0', energies)

    grad = casscf.compute_gradient(2, traj_id=traj_id)
    print(f'Trajectory {traj_id} gradient root 2 at step 0', grad)

for traj_id, geom in enumerate(geom_step1):
    # Second geometry reuses only this trajectory's previous CASSCF/HF state.
    casscf.set_geom_and_run_hf(geom, traj_id=traj_id)

    energies = [casscf.compute_energy(root, traj_id=traj_id) for root in range(NSTATES)]
    print(f'Trajectory {traj_id} energies at step 1', energies)

    grad = casscf.compute_gradient(2, traj_id=traj_id)
    print(f'Trajectory {traj_id} gradient root 2 at step 1', grad)

    overlap = casscf.time_overlap_matrix(NSTATES, traj_id=traj_id)
    print(f'Trajectory {traj_id} time-overlap matrix', overlap)

assert isinstance(casscf, ElectronicStructureStrategy)
