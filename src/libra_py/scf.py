#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  This module contains a number of functions to  compute properties of interest
#  from the converged electronic structure calculations
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import data_outs
from . import pdos

def takeSecond(elem):
    return elem[1]

def spectrum(ham, T_file = "T_mo.dat", spec_file = "spectrum.txt"):
    """
    ham - is a listHamiltonian object, that contains basis and converged MO info
    T_file - the name of the file to store the transition dipole moment
    spec_file - the name of the file to store the energy - oscillator strength pairs (sorted)

    Compute the electronic absorption spectrum, based on the transition dipole moments
    
    """

    bas = ham.basis_ao
    es = ham.get_electronic_structure()

    homo = es.Nocc_alp - 1

    C = es.get_C_alp()
    E = es.get_E_alp()

    nao, nmo = C.num_of_rows, C.num_of_cols

    Tao_x = MATRIX(nao, nao)
    Tao_y = MATRIX(nao, nao)
    Tao_z = MATRIX(nao, nao)

    for a in range(0,nao):
        for b in range(0,nao):
            mu = transition_dipole_moment(bas[a], bas[b])
            Tao_x.set(a,b, mu.x)
            Tao_y.set(a,b, mu.y)
            Tao_z.set(a,b, mu.z)

    Tmo_x = C.T() * Tao_x * C
    Tmo_y = C.T() * Tao_y * C
    Tmo_z = C.T() * Tao_z * C

    Tmo = MATRIX(nmo, nmo)

    res = []
    for i in range(0,homo+1):   # All unoccupied orbitals
        for j in range(i+1, nmo):   #  All occupied orbitals
            dE = E.get(j,j)-E.get(i,i)
            X = Tmo_x.get(i,j)
            Y = Tmo_y.get(i,j)
            Z = Tmo_z.get(i,j)
            Tmo.add(i,j, (2.0/3.0)*dE*( X*X + Y*Y + Z*Z) )
            res.append([Tmo.get(i,j), dE])


    # sort list with key
    sorted_res = sorted(res, key=takeSecond)

    data_outs.show_matrix_splot(Tmo, T_file)        

    f = open(spec_file,"w") 
    for it in sorted_res:
        f.write("%8.5f %8.5f\n" % (it[1], it[0]))
    f.close()
  


def print_orbitals(ham, syst, orbs, prefix = "", grid = [40, 40, 40]):
    """
    This function will print out the wavefunctions of the converged SCF calculations

    ham - the listHamiltonian object containing the basis and MO info
    syst - chemical system object, containing the coordinates info
    obrs - a list of orbitals to be processes (list of integers)
    prefix - a prefix of the filenames where the wavefunctions be printed
    grid - 3 integers defining the coarsening level of the cube file

    """

    prms = Control_Parameters()
    prms.nx_grid, prms.ny_grid, prms.nz_grid = grid[0], grid[1], grid[2]
    prms.charge_density_prefix = prefix+"_orbital_"
    prms.orbs = Py2Cpp_int(orbs) 

    es = ham.get_electronic_structure()
    charge_density( es, syst, ham.basis_ao, prms)



def print_pdos(ham, syst, projections, prefix = "pdos", emin = -35.0, emax = 35.0, de = 0.1, outfile="pdos.txt"):
    """
    ham - listHamiltonian object that contains info about atom2ao mapping, AO basis, and MOs
    syst - contains the atomic coordinates
    projections - which types of DOS projections to compute. Format like in [["tot",range(0,syst.Number_of_atoms)]]
    prefix - the directory name (to be created if not present) where all the files will go

    """

    prms = Control_Parameters()
    prms.dos_prefix = prefix+"/"

    os.system("mkdir " + prefix)

    es = ham.get_electronic_structure()
    compute_dos( es, syst, ham.basis_ao, prms, ham.atom_to_ao_map)
    proj_dos.pdos(emin, emax, de, projections, prefix+"/_alpha_wfc_atom", outfile, es.Nelec)





