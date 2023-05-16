#*********************************************************************************
#* Copyright (C) 2018-20203 Mahsa Mofidi, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: LAMMPS_methods
   :platform: Unix, Windows
   :synopsis: This module implements various functions for computing with LAMMPS
.. moduleauthor:: Alexey V. Akimov

"""

import os
import math
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def compute_dynmat(lmp, filename, atoms, dr, opt=1):
    """
    The function to return the dynamic matrix for a system using LAMMPS

    Note: 
        To use this code, one needs to have pylammps installed

    Args:
        lmp ( lammps.lammps() object ): an instance of PyLAMMPS - you need to create if first
        filename ( string ): LAMMPS input filename defining the system and the interactions
            (importantly! assuming everything is in ATOMIC UNITS - "units electron"). Furthermore,
            this file should define the minimization procedure e.g. "minimize  0.0 1.0e-8 1000 100000"
        atoms ( list of integers ): indices of the atoms for which to compute Hessian
        dr ( double ): the magnitude of the displacement [ units: Bohr ]

    Returns: 
        tuple: ( D, H, R, M, E ), where:

            * D ( MATRIX(ndof, ndof) ): a sub-block of a dynamic matrix (Hessian/sqrt(m_ij)) 
            * H ( MATRIX(ndof, ndof) ): a sub-block of a Hessian matrix [ units: a.u.]
            * r ( MATRIX(ndof, 1) ): coordinates of all atoms [ units: Bohr ]
            * M ( MATRIX(ndof, 1) ): masses of all atoms [ units: a.u. ]
            * E ( list of ndof ints ): types of atoms - later can be mapped to atom names
            * opt ( int ): flag to use all atoms or only the pupplied sub-set

                - opt = 0: use all atoms
                - opt = 1: use only the supplied subset [ default ]

        Example:
            >> import lammps 
            >> lmp = lammps.lammps()
            >> D, H, R, M, E = compute_dynmat(lmp, "in.lammps", [0,1,2], 1e-3)
            >> print "Atom types are:", E
            >> print "Masses are:"; M.show_matrix()
            >> print "Coordinates are:"; R.show_matrix()
            >> print "Hessian matrix is:"; H.show_matrix()
            >> print "Dynamical matrix is:"; D.show_matrix()
    """

    # run infile all at once
    lmp.file(filename)
   
    nlocal = lmp.extract_global("nlocal", 0)
    r = lmp.numpy.extract_atom_darray("x", nlocal, dim=3)
    F = lmp.numpy.extract_atom_darray("f", nlocal, dim=3)

    ntypes = lmp.extract_global("ntypes", 0)
    mass = lmp.numpy.extract_atom_darray("mass", ntypes+1)
    atype = lmp.numpy.extract_atom_iarray("type", nlocal)

    if opt==0:
        atoms = range(nlocal)

    nat = len(atoms)
    ndof = 3 * nat
    H = MATRIX(ndof, ndof)
    D = MATRIX(ndof, ndof)
    R = MATRIX(ndof, 1)
    M = MATRIX(ndof, 1)
    E = []


    # Coordinates
    for at in atoms:
        a = atoms.index(at)
        for pa in [0,1,2]:
            i = 3*a + pa
            R.set(i,0, r[at][pa])

    # Masses and types
    for at in atoms:
        a = atoms.index(at)
        E.append(atype[at][0])
        for pa in [0,1,2]:
            i = 3*a + pa      
            M.set(i,0, mass[atype[at][0]][0])


    # Hessian and Dynamtrix
    for at in atoms:    
        a = atoms.index(at)
        for pa in [0,1,2]:
            i = 3*a + pa 

            # Theory:
            # dH_ij = (F_i(R+R_j) - F_i(R-R_j)) / 2*dR

            # Displace the atom a by dr in direction pa
            # so we evaluate all forces at R+dR[at][pa] 
            # update the forces
            r[at][pa] = r[at][pa] + dr
            lmp.command("run 0")        

            for bt in atoms:
                b = atoms.index(bt)
                for pb in [0,1,2]:
                    j = 3*b+pb
                    H.add(j,i, F[bt][pb])

            # Displace the atom a by 2*dr in direction -pa
            # so we evaluate all forces at R-dR[at][pa]
            # update the forces
            r[at][pa] = r[at][pa] - 2.0*dr       
            lmp.command("run 0")
        
            for bt in atoms:
                b = atoms.index(bt)
                for pb in [0,1,2]:
                    j = 3*b+pb
                    H.add(j,i, -1.0*F[bt][pb])

            # Return the atom a to its original position
            r[at][pa] = r[at][pa] + dr

    H = (0.5/dr)*H

    for i in range(0,ndof):
        for j in range(0,ndof):
            D.set(i,j, H.get(i,j)/math.sqrt(M.get(i)*M.get(j)))

    return D, H, R, M, E



"""

"""






