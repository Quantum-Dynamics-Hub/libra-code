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
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def load_mol():
    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")

    # Create molecular system and initialize the properties
    syst0 = System();  LoadMolecule.Load_Molecule(U, syst0, "system.pdb", "true_pdb2")
    syst0.determine_functional_groups(0)  # do not assign rings

    # Define the groups: note - the indexing of atoms and fragments starts from 1, not 0
    syst0.GROUP_ATOMS(range(1,61), 1)    # C60
    syst0.GROUP_ATOMS(range(61,241), 2)  # graphene 

    syst0.init_fragments()
    print "Number of atoms in the system = ", syst0.Number_of_atoms
    print "Number of bonds in the system = ", syst0.Number_of_bonds
    print "Number of angles in the system = ", syst0.Number_of_angles
    print "Number of dihedrals in the system = ", syst0.Number_of_dihedrals
    print "Number of impropers in the system = ", syst0.Number_of_impropers
    print "Number of fragments in the system = ", syst0.Number_of_fragments

    return syst0



def test_rotations(case):

    syst0 = load_mol()
     
    u, pivot_pt = None, None    

    if case==100:
        # Rotation around the center of mass around x axis
        pivot_pt = syst0.Fragments[0].Group_RB.rb_cm
        u = VECTOR(1.0, 0.0, 0.0)

    elif case==101:  
        # Rotation around the center of mass around y axis
        pivot_pt = syst0.Fragments[0].Group_RB.rb_cm
        u = VECTOR(0.0, 1.0, 0.0)

    elif case==102:  
        # Rotation around the center of mass around z axis
        pivot_pt = syst0.Fragments[0].Group_RB.rb_cm
        u = VECTOR(0.0, 0.0, 1.0)



    elif case==200:  
        # Rotation around an external center around x axis
        pivot_pt = VECTOR(0.0, 0.0, 0.0)
        u = VECTOR(1.0, 0.0, 0.0)

    elif case==201:  
        # Rotation around an external center around y axis
        pivot_pt = VECTOR(0.0, 0.0, 0.0)
        u = VECTOR(0.0, 1.0, 0.0)

    elif case==202:  
        # Rotation around an external center around y axis
        pivot_pt = VECTOR(0.0, 0.0, 0.0)
        u = VECTOR(0.0, 0.0, 1.0)


    elif case==300:  
        # Rotation around an atom around x axis
        pivot_pt = syst0.Atoms[0].Atom_RB.rb_cm
        u = VECTOR(1.0, 0.0, 0.0)

    elif case==301:  
        # Rotation around an atom around y axis
        pivot_pt = syst0.Atoms[0].Atom_RB.rb_cm
        u = VECTOR(0.0, 1.0, 0.0)

    elif case==302:  
        # Rotation around an atom around z axis
        pivot_pt = syst0.Atoms[0].Atom_RB.rb_cm
        u = VECTOR(0.0, 0.0, 1.0)


    elif case==400:  
        # Rotation around an atom around the bond between two atoms
        pivot_pt = syst0.Atoms[0].Atom_RB.rb_cm
        u = syst0.Atoms[0].Atom_RB.rb_cm - syst0.Atoms[1].Atom_RB.rb_cm


   
    f = open("rot-case-%i.xyz" % (case),"w");  f.close()  

    dphi = 1.0
    frag_id = 1

    for i in xrange(360): 
        syst0.ROTATE_FRAGMENT(dphi, u, frag_id,  pivot_pt)
        syst0.print_xyz("rot-case-%i.xyz" % (case), i)



def test_translations(case):

    syst0 = load_mol()
     
    u = None

    if case==100:
        u = VECTOR(1.0, 0.0, 0.0)

    elif case==101:
        u = VECTOR(0.0, 1.0, 0.0)

    elif case==102:
        u = VECTOR(0.0, 0.0, 1.0)



    f = open("tr-case-%i.xyz" % (case),"w");  f.close()  

    dr = 0.1
    frag_id = 1

    for i in xrange(100): 
        syst0.TRANSLATE_FRAGMENT(dr, u, frag_id)  
        syst0.print_xyz("tr-case-%i.xyz" % (case), i)




def test_mixed(case):

    syst0 = load_mol()
     
    u, pivot_pt, tdir = None, None, None    

    if case==100:
        u = VECTOR(0.0, 1.0, 0.0)
        tdir = VECTOR(1.0, 0.0, 0.0)



    f = open("mixed-case-%i.xyz" % (case),"w");  f.close()  

    dphi = 3.6
    dr = 0.1
    frag_id = 1

    for i in xrange(100): 
        syst0.TRANSLATE_FRAGMENT(dr, tdir, frag_id)  
        pivot_pt = syst0.Fragments[0].Group_RB.rb_cm
        syst0.ROTATE_FRAGMENT(dphi, u, frag_id,  pivot_pt)

        syst0.print_xyz("mixed-case-%i.xyz" % (case), i)




#================= Main part =================
for case in [100, 101, 102, 200, 201, 202, 300, 301, 302, 400]:
    test_rotations(case)
 
for case in [100, 101, 102]:
    test_translations(case)

for case in [100]:
    test_mixed(case)

