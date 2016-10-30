#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import math
import sys
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


# Conversion parameters
Angst_to_Bohr = 1.889725989  

#define constants
a	= 4.07825
Nx	= 6
Ny	= 6
Nz	= 2

#define lattice vector (unit cell)
t1 = VECTOR(0.5*a*math.sqrt(2.0),        0.0,                         0.0 )
t2 = VECTOR(0.5*a*math.sqrt(1.0/2.0),    0.5*a*math.sqrt(3.0/2.0),    0.0 )
t3 = VECTOR(0.0,                         0.0,                         a*math.sqrt(3.0) )

# and the transformation matrix
O = MATRIX(3,3)
O.set(0,0,t1.x); O.set(0,1,t2.x); O.set(0,2,t3.x)
O.set(1,0,t1.y); O.set(1,1,t2.y); O.set(1,2,t3.y)
O.set(2,0,t1.z); O.set(2,1,t2.z); O.set(2,2,t3.z)

#define atomic basis
Au = []
Au.append( VECTOR( 0.0,      0.0,       0.0   ) )
Au.append( VECTOR(1.0/3.0,  1.0/3.0,   1.0/3.0) )
Au.append( VECTOR(2.0/3.0,  2.0/3.0,   2.0/3.0) )



#Compute the number of atoms,using input information: Nx, Ny, Nz
Nat = len(Au)*Nx*Ny*Nz   

# Create a Universe and the empty system
U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")
syst = System()


#compute atomic coordinates and print them into the file
myFile = open('atoms.xyz','w')
myFile.write('%5i\n' % Nat)
myFile.write('\n')

for nx in range(Nx):
   for ny in range(Ny):
      for nz in range(Nz):

         for i in range(0,3):
            # Crystal coordinates = unit cell atomic basis + integer translations
            au_i_cryst = VECTOR( Au[i].x + nx, Au[i].y + ny, Au[i].z + nz )

            #transform crystal coordinates to Cartesian
            au_i_cart = O * au_i_cryst
            
            # Print coordinates
            myFile.write('Au\t% 16.5f\t% 16.5f\t% 16.5f\n' % (au_i_cart.x, au_i_cart.y, au_i_cart.z))  # in Angstrom


            # Create Au atoms and add them to the chemical system
            at = Atom(U, {"Atom_element": "Au", 
                          "Atom_cm_x": au_i_cart.x * Angst_to_Bohr,
                          "Atom_cm_y": au_i_cart.y * Angst_to_Bohr, 
                          "Atom_cm_z": au_i_cart.z * Angst_to_Bohr
                         } )
            syst.CREATE_ATOM(at)  # internally, we assume Bohrs, so the conversion is made up there

#close the output file
myFile.close()


# Compute the supercell translation vectors
T1 = Nx*t1; T2 = Ny*t2; T3 = Nz*t3

# Convert lattice vectors to Bohrs
T1 *= Angst_to_Bohr
T2 *= Angst_to_Bohr
T3 *= Angst_to_Bohr

print "Lattice vectors"
print "T1 = ", T1.x, T1.y, T1.z
print "T2 = ", T2.x, T2.y, T2.z
print "T3 = ", T3.x, T3.y, T3.z

syst.init_box(T1,T2,T3);
syst.print_ent("au.pdb");




