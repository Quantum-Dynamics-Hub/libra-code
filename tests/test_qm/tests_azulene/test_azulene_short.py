#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
# Tutorial: SCF computations are hidden - use built-in function
###################################################################
import os
import math
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


#=========== STEP 1:  Create Universe and populate it ================
U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.dat", 1)

#=========== STEP 2:  Create system and load a molecule ================
syst = System()
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/azulene.pdb", "iqmol_pdb")
#Load_Molecule(U, syst, os.getcwd()+"/benzene.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/c.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/c2.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/ch4.pdb", "pdb_1")

print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)

#=========== STEP 3: Compute electronic structure ================
lstHamQM = listHamiltonian_QM("control_parameters_eht.dat", syst)
lstHamQM.compute_scf(syst)
el_str = lstHamQM.get_electronic_structure()

print "Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
for j in xrange(el_str.Norb):
    print "%12.8f   %12.8f  %12.8f   %12.8f" %(el_str.get_bands_alp(j), el_str.get_occ_alp(j), el_str.get_bands_bet(j), el_str.get_occ_bet(j) )

#=========== STEP 4: Compute charge density for HOMO and LUMO ================
# Compute homo index
homo = el_str.Nocc_alp - 1   # index of the HOMO orbital

prms = Control_Parameters()
prms.orbs = Py2Cpp_int([homo, homo+1]) 
prms.nx_grid, prms.ny_grid, prms.nz_grid = 40, 40, 40
prms.charge_density_prefix = "char_dens/"

charge_density( el_str, syst, lstHamQM.basis_ao, prms)

