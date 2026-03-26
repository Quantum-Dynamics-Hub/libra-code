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

from pyscf import gto
from pyscf2.implementations.casscf import CASSCF
from pyscf2.interfaces import ElectronicStructureStrategy
import numpy as np

from pyscf2.interfaces import MolecularGeometry

geom1 = MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0,0.0,0.0],[0.0,0.0,0.7746]]))
geom2 = MolecularGeometry(atom_labels=['He', 'H'], coords_angstrom=np.array([[0.0,0.0,0.0],[0.0,0.0,1.0]]))



# instantiate the CASSCF strategy class with parameters for the test
casscf = CASSCF(norbcas=2, nelecas=2, nroots=3, basis='sto-3g', charge=1)

# wrap geometry set and HF in set_geom_and_run_hf method
casscf.set_geom_and_run_hf(geom1)

# 2) energy for root 0,1,2
energies1 = [casscf.compute_energy(root) for root in range(3)]
print('Energies at geom1', energies1)

# 3) gradient for root 2
grad1 = casscf.compute_gradient(2)
print('Gradient root 2 geom1', grad1)

# 4) set mol2 geometry and run HF
casscf.set_geom_and_run_hf(geom2)

# 5) energy for root 0,1,2
energies2 = [casscf.compute_energy(root) for root in range(3)]
print('Energies at geom2', energies2)

# 6) gradient for root 2
grad2 = casscf.compute_gradient(2)
print('Gradient root 2 geom2', grad2)

# 7) time-overlap matrix for 3 roots
overlap = casscf.time_overlap_matrix(3)
print('Time-overlap matrix', overlap)

assert isinstance(casscf, ElectronicStructureStrategy)


