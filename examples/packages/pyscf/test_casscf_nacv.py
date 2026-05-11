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
.. module:: pyscf.implementations.test_lif_scan
   :platform: Unix, Windows
   :synopsis: Test for LiF PES scan using universal ES interface.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

import sys
from pathlib import Path
import numpy as np

#def _prepend_repo_root() -> None:
#    file_path = Path(__file__).resolve()
#    for parent in file_path.parents:
#        if (parent / "interface" / "__init__.py").is_file():
#            repo_root = str(parent)
#            if repo_root not in sys.path:
#                sys.path.insert(0, repo_root)
#            return
#
#if __name__ == "__main__" and __package__ is None:
#    _prepend_repo_root()

#try:
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry

#except ModuleNotFoundError as exc:
#    if exc.name not in {"interface", "interface.implementations.casscf", "interface.interfaces"}:
#        raise
#    from libra_py.packages.pyscf.implementations.casscf import CASSCF
#    from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry

NSTATES = 2
DISTANCE_START_BOHR = 7.0
DISTANCE_STOP_BOHR = 15.0
DISTANCE_STEP_BOHR = 1.0

basis_dict = {'Li': 'sto-3g', 'F': '6-311+g*'}
cas_list = [4, 7, 11, 14, 17]# 0-indexed: 3 (F2pz), 6 (Li2s), 10 (F5pz), 13 (F5s), 16 (F4pz)

distances = np.arange(DISTANCE_START_BOHR, DISTANCE_STOP_BOHR, DISTANCE_STEP_BOHR, dtype=np.float64)

# Initialize the CASSCF strategy
casscf = CASSCF(
    norbcas=5, 
    nelecas=2, 
    nroots=NSTATES, 
    basis=basis_dict, 
    unit='Bohr', 
    charge=0,
    cas_list=cas_list
)

for step, d in enumerate(distances):
    print(f"\n========== Distance: {d:.2f} Bohr ==========")
    # Though the parameter is named 'coords_angstrom', passing raw units matching
    # the target unit (Bohr) works properly because we set unit='Bohr' in CASSCF.
    geom = MolecularGeometry(
        atom_labels=['Li', 'F'], 
        coords_angstrom=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, d]])
    )
    
    casscf.set_geom_and_run_hf(geom)
    
    energies = [casscf.compute_energy(root) for root in range(NSTATES)]
    print(f"Energies: {energies}")
    
    # Calculate pairwise NACVs (requires use_etfs=False to match your previous scan logic)
    nacv = casscf.compute_nac_vectors(use_etfs=False)
    
    # Print F z-direction NACV matrix
    F_ATOM_INDEX = 1
    Z_AXIS_INDEX = 2
    f_z_nacv = nacv[:, :, F_ATOM_INDEX, Z_AXIS_INDEX]
    print(f"F z-direction NACV matrix:\n{f_z_nacv}\n")
    
    if step > 0:
        overlap = casscf.time_overlap_matrix(NSTATES)
        print("Time-overlap matrix with previous geometry:")
        print(overlap)
