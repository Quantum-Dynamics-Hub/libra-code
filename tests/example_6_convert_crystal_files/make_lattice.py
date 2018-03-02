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
"""
  Prepare a perovskite supercell from sample input files

"""


import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


# Conversion from Angstrom to Bohr units, so that you can input everything 
# in the convenient Angstrom units
Angst = 1.889725989     

Nx, Ny, Nz = 1, 1, 1

#prefix = "perov_majet"
prefix = "perov_weili"

# Read the minimal unit cell information from the file
L0 = []
R0 = []

if prefix=="perov_majet":

    cell_size = 16.9843 / Angst   #  cubit unit cell size, in Angstrom

    a = VECTOR(16.9843, 0.0, 0.0) # convert to a.u.
    b = VECTOR(0.0, 16.9843, 0.0) # convert to a.u.
    c = VECTOR(0.0, 0.0, 16.9843) # convert to a.u.

    L0, R0raw = build.read_xyz_crystal(prefix+".xyz", a, b, c)  # in Bohr

    for r in R0raw:
        R0.append( r )   # coordinates in a.u.
        print r.x, r.y, r.z



elif prefix=="perov_weili":

    L0, R0raw = build.read_xyz(prefix+".xyz")  # in Bohr
    Nx, Ny, Nz = 2, 1, 1

    a = VECTOR(8.8009, 0.0,    0.0) * Angst # convert to a.u.
    b = VECTOR(0.0,    8.009,  0.0) * Angst # convert to a.u.
    c = VECTOR(0.0,    0.0,   12.6857) * Angst # convert to a.u.

    for r in R0raw:
        R0.append( r )   # coordinates in a.u.
        print r.x, r.y, r.z



L1, R1 = build.generate_replicas_xyz2(L0, R0, a, b, c, Nx, Ny, Nz)

masses = []
for i in xrange(len(L1)):
    masses.append(1.0)


# Create an empty chemical system object
syst = System()
build.add_atoms_to_system(syst, L1, R1, VECTOR(1.0, 1.0, 1.0), masses, "elements.dat")
syst.init_box(a*Nx,b*Ny,c*Nz)

# And print it out in the .ent format. Note, this will only print out the coordinates, and not the connectivities
syst.print_ent(prefix+".ent") 


