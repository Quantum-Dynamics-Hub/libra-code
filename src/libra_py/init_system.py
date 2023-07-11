#*********************************************************************************
#* Copyright (C) 2017-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: init_system
   :platform: Unix, Windows
   :synopsis: This module implements functions to initialize System() objects 
.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import LoadPT  # for Load_PT


def init_system(label, R, g, rnd, T, sigma, df, data_file="elements.dat"):
    """This function creates and instance of the System object which will be used in classical MD.

    Args:
        label ( list of Nat strings ): Atom names
        R ( list of Nat VECTORs ): Coordinates of all atoms [ units: Bohr ]
        g ( list of Nat VECTORs ): Gradients on all atoms [ units: Ha/Bohr ]
        rnd ( Random ): Random number generator object
        T ( double ): Target temperature used to initialize momenta of atoms [ units: K ]
        sigma ( double ): The magnitude of a random displacement of each atom from its center [ units: Bohr ]
        df ( double ): Controls additional (debug) printing
        data_file ( string ): The file containing information about elements.
 
    Returns:
        ( System ): syst : The System object that contains all the geometry information
 
    """

    # Create Universe and populate it
    U = Universe();   LoadPT.Load_PT(U, data_file, 0)
    syst = System()

    sz = len(label)
    for i in range(0,sz):
        atom_dict = {"Atom_element": label[i]}

        if df:
            print("CREATE_ATOM ",atom_dict["Atom_element"])

        at = Atom(U, atom_dict)
        r = VECTOR( R[i].x + sigma*rnd.normal(), R[i].y + sigma*rnd.normal(), R[i].z + sigma*rnd.normal() )
        print(r, r.x, r.y, r.z)
        at.Atom_RB.set_position( r )
        at.Atom_RB.set_force( VECTOR(-g[i].x, -g[i].y, -g[i].z) )
        at.is_Atom_RB = 1
        syst.CREATE_ATOM(at)

    if df:
        print("Number of atoms in the system = ", syst.Number_of_atoms)
        syst.show_atoms()

    # Initialize random velocity at T(K) using normal distribution
    syst.init_atom_velocities(T,rnd)

    return syst

