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
  This file shows how to build a CdSe quantum dot (QD) of a given size.
  Here, we don't cap the surface atoms, but still
  generate the connectivity map for the final QD to save the file in the .ent format
  for futher use in MD simulations

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


# Unit cell parameters
# Periodic translations of the minimal unit cell (provided as an input file)
a = VECTOR(5.833, 0.0, 0.0) * Angst
b = VECTOR(0.0, 5.833, 0.0) * Angst
c = VECTOR(0.0, 0.0, 5.833) * Angst

# Read the minimal unit cell information from the file
L0, R0 = build.read_xyz("CdSe.xyz")


# Now, lets create a super-cell by replicating the minimal cell in 3 directions
# We need large enough super-cell so that the QDs we are going to create fits inside it
# so, replicate 6 times in each direction
Nx, Ny, Nz = 6, 6, 6 

# This will generate the labels and coordinates of the new super-system
L1, R1 = build.generate_replicas_xyz2(L0, R0, a, b, c, Nx, Ny, Nz)


# Next, we'll be connecting atoms in the super-cell to each other, based on their separation distances
# and the maximal coordination numbers accepted
# At this point, we are going to define the maximal coordination number of each atom type present in the 
# system (here, Cd and Se), as well as for each additional elements - e.g. those which we are going to cap
# the QD with (here, H)

PT_coord = { "Cd": 4, "Se": 4, "H":1 }

MaxCoord1 = []
masses = []
for i in L1:
    MaxCoord1.append(PT_coord[i])    # maximal coordination number of each atom in super-cell
    masses.append(1.0)               # masses of each atom in the super-cell, do not matter here



# The parameter that defines the atomic connectivities: if the atoms are separated
# the the distance shorter than this, and if their actual coordination number has not
# rached the maximal one yet, we will assume they are covelently bonded (in FF there will
# be a harmonic potential term for instance)
Rbond = 3.001 * Angst


# The radius of the QD 
Rdot = 7.0 * Angst

# Finally, we can carve a QD out of the super-cell. All the atoms that are outside the shepre of Rdot radius
# are removed. The L and R variables will contain the labels and coordinates of the final (QD) atoms.
L, R = build.crop_sphere_xyz2(L1, R1, Rdot)


# Now, we are going to auto-connect all the atoms in the final QD
# Define the maximal coordination for the atoms 
MaxCoord2 = []
for i in L:
    MaxCoord2.append(PT_coord[i])   # maximal coordination number of each atom in the resulting QD

# Determine the connectivities of the new set of atoms
res, line, pairs = autoconnect.autoconnect_pbc(R, MaxCoord2, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0)
print "Autodetermined connections: "
print line


# Create an empty chemical system object
syst = System()

# Add all the atoms of the final QD to the chemical system object
build.add_atoms_to_system(syst, L, R, VECTOR(1.0, 1.0, 1.0), masses, "elements.dat")

# And print it out in the .ent format. Note, this will only print out the coordinates, and not the connectivities
syst.print_ent("CdSe-qd.ent") 

# We print out the connectivities afterwards manually
f = open("CdSe-qd.ent", "a")
f.write("END\n")
f.write(line)
f.write("END\n")
f.close()

