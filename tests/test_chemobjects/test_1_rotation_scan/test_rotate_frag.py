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
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def main():
    rnd = Random()
    #--------------------- Initialization ----------------------

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")

    # Create molecular system and initialize the properties
    syst0 = System();  LoadMolecule.Load_Molecule(U, syst0, "c2h4.ent", "pdb")


    syst0.determine_functional_groups(0)  # do not assign rings
    syst0.init_fragments()
    print "Number of atoms in the system = ", syst0.Number_of_atoms
    print "Number of bonds in the system = ", syst0.Number_of_bonds
    print "Number of angles in the system = ", syst0.Number_of_angles
    print "Number of dihedrals in the system = ", syst0.Number_of_dihedrals
    print "Number of impropers in the system = ", syst0.Number_of_impropers


    # Generate configurations
    print "Printing coordinates of all atoms" 
    a = VECTORList()
    for i in xrange(syst0.Number_of_atoms):
        a.append( VECTOR( syst0.Atoms[i].Atom_RB.rb_cm.x, syst0.Atoms[i].Atom_RB.rb_cm.y, syst0.Atoms[i].Atom_RB.rb_cm.z ) )
        print i, a[i]

    v1 = a[3] - a[0]  # C-C direction
    v2 = a[0] - a[2]
    
    v3 = VECTOR(v2.x, 0, 0) 
    v4 = VECTOR(0, v2.y, 0)

    f = open("c2h4_0.xyz","w");  f.close()
    fr_num = 1  # ID of the fragment we are going to rotate
    for i in xrange(10):
 
        syst0.ROTATE_FRAGMENT(18.0, v3, fr_num, 0)                

        syst0.print_xyz("c2h4_0.xyz",i)

 
main()
