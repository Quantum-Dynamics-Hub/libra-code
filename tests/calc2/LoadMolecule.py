#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import re
import os
import sys
import math


# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/chemobjects")

from cygmmath import *
from cygchemobjects import *

Angst_to_Bohr = 1.889725989  


def Load_Molecule(univ,syst,mol_file,format):
#-------- Load molecular system -------------
# Specify format manually (this gives flexibility)
# of format recognizable by program
# In this particular implementation my .ent format is
# loading


######### Assume coordinates are given in Angstroms #############

#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'    
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '"((\w|\W)+)"'
    compINT = re.compile(INT)


#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_keyword = '(?P<Atom_keyword>'+'HETATM'+')'+SP
    pAtom_id      = '(?P<Atom_id>'+DOUBLE+')'+SP    
    pAtom_element = '(?P<Atom_element>'+WORD+')'+SP
    pAtom_chain   = '(?P<Atom_chain>'+WORD+')'+SP
    pAtom_id1     = '(?P<Atom_id1>'+DOUBLE+')'+SP
    pAtom_x_coord = '(?P<Atom_x_coord>'+DOUBLE+')'+SP
    pAtom_y_coord = '(?P<Atom_y_coord>'+DOUBLE+')'+SP
    pAtom_z_coord = '(?P<Atom_z_coord>'+DOUBLE+')'+SP
    pAtom_type    = '(?P<Atom_type>'+INT+')'+SP
    pAtom_charge  = '(?P<Atom_charge>'+DOUBLE+')'+SP
    pAtom_C60     = '(?P<Atom_C60>'+NINT+')'+SP

    if format=="pdb":
        Atom_Record = pAtom_keyword + pAtom_id + pAtom_element + pAtom_id1 + pAtom_x_coord + pAtom_y_coord + pAtom_z_coord + pAtom_charge 
    elif format=="pdb_1":
        Atom_Record = pAtom_keyword + pAtom_id + pAtom_element + pAtom_id1 + pAtom_x_coord + pAtom_y_coord + pAtom_z_coord + pAtom_chain
    elif format=="true_pdb":
        Atom_Record = pAtom_keyword + pAtom_id + pAtom_element + pAtom_mol + pAtom_chain + pAtom_id1 + pAtom_x_coord + pAtom_y_coord + pAtom_z_coord + pAtom_occ + pAtom_charge 



    pBond_keyword = '(?P<Bond_keyword>'+'CONECT'+')'+SP

    Bond_Record = pBond_keyword + '('+pAtom_id+')+'

    
    pFrag_keyword1 = '(?P<Frag_keyword1>'+'GROUP'+')'+SP
    pFrag_id       = '(?P<Frag_id>'+INT+')'+SP
    pFrag_keyword2 = '(?P<Frag_keyword2>'+'FRAGNAME'+')'+SP
    pFrag_fragname = '(?P<Frag_fragname>'+ID+')'+SP
    pFrag_atom1name= '(?P<Frag_atom1name>'+ID+')'+SP
    pFrag_atom2name= '(?P<Frag_atom2name>'+ID+')'+SP
    pFrag_atom3name= '(?P<Frag_atom3name>'+ID+')'+SP

    Fragment_Record = pFrag_keyword1 + '('+pFrag_id+')+' + \
                      pFrag_keyword2 + pFrag_fragname + \
                      pFrag_atom1name + pFrag_atom2name + \
                      pFrag_atom3name
  

    f   = open(mol_file,'r')
    A = f.readlines()
    D = C = B = A
    f.close()

    #---------- Create atoms ----------------
    for a in A:        
        m1 = re.search(Atom_Record,a)
        if m1!=None:           
            atom_dict = {}

            atom_dict["Atom_element"] = a[m1.start('Atom_element'):m1.end('Atom_element')]
            atom_dict["Atom_cm_x"] = float(a[m1.start('Atom_x_coord'):m1.end('Atom_x_coord')]) * Angst_to_Bohr
            atom_dict["Atom_cm_y"] = float(a[m1.start('Atom_y_coord'):m1.end('Atom_y_coord')]) * Angst_to_Bohr
            atom_dict["Atom_cm_z"] = float(a[m1.start('Atom_z_coord'):m1.end('Atom_z_coord')]) * Angst_to_Bohr

            if format=="pdb" or format=="true_pdb":
                atom_dict["Atom_charge"] = float(a[m1.start('Atom_charge'):m1.end('Atom_charge')])
#            atm.Atom_ff_int_type     = int(float(a[m1.start('Atom_type'):m1.end('Atom_type')]))
#            atm.Atom_is_C60_CT       = int(a[m1.start('Atom_C60'):m1.end('Atom_C60')])
 
            print "CREATE_ATOM ",atom_dict["Atom_element"]
            syst.CREATE_ATOM( Atom(univ,atom_dict)  ) 
    
    #--------- Create bonds ------------------
    print len(B)
    for b in B:
        m2 = re.search(Bond_Record,b)
        lst = []
        if m2!=None:
            lst  = compINT.findall(m2.group(0))
#            if (int(float(lst[0][0])) == 10):
#               exit(0)
        i = 1   
        while i < len(lst):
            print 'LINK_ATOMS ',lst[0][0], ' and ', lst[i][0],'...'
            syst.LINK_ATOMS(int(float(lst[0][0])),int(float(lst[i][0])))
            i = i + 1
        print "----------------------------"

    #----------- Group atoms --------------------
    j = 1
    for fr in C:
        m3 = re.search(Fragment_Record,fr)        
        lst = []
        if m3!=None:
            print m3.group()
            lst = compINT.findall(m3.group())

            gr_atoms = []
            i = 0
            while i<len(lst):
                gr_atoms.append(int(float(lst[i][0])))
                i = i + 1
#            if len(gr_atoms)>0:
            print "GROUP_ATOMS",gr_atoms," to form fragment with Group_id = ",j
            syst.GROUP_ATOMS(gr_atoms,j)
#            syst.show_fragments()
            j = j + 1

    #----------------------------------------------------------------------

