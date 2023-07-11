#*********************************************************************************
#* Copyright (C) 2015-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: LoadMolecule
   :platform: Unix, Windows
   :synopsis: This module implements functions for loading data into a Chemobject 
       objects - by reading formatted data files. 
.. moduleauthor:: Alexey V. Akimov

"""
import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from .regexlib import *
from . import units


def Load_Molecule(univ,syst,mol_file,format, verbosity=0):
    """Load molecular system from various formats
 
    Specify format manually (this gives flexibility) of format recognizable by program

    Args:
        univ ( Universe object ): Contains the basic information about chemical elements
        syst ( System ): The Chemobjects object that represents molecular system - this is what we construct
        mol_file ( string ): The name of the file containing the molecular structure
        format ( string ): The name of the format according to which the file "mol_file" is assumed to be formatted
            Available options are:
 
            * pdb
            * pdb_1
            * true_pdb
            * true_pdb2
            * xyz
            * iqmol_pdb

    Returns:
        None: but the `syst` object is modified to add the atoms and bonds, and to group atoms together

    Note:
        We assume that the coordinates in the files read are given in Angsrom.

        However, the syst stores this data in the atomic units (Borh), so conversion happens

    """

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_keyword = '(?P<Atom_keyword>'+'HETATM'+')'+SP
    pAtom_keyword1= '(?P<Atom_keyword>'+'ATOM'+')'+SP
    pAtom_id      = '(?P<Atom_id>'+DOUBLE+')'+SP    
    pAtom_element = '(?P<Atom_element>'+ID+')'+SP
    pAtom_mol     = '(?P<Atom_mol>'+WORD+')'+SP
    pAtom_chain   = '(?P<Atom_chain>'+WORD+')'+SP
    pAtom_id1     = '(?P<Atom_id1>'+DOUBLE+')'+SP
    pAtom_type    = '(?P<Atom_type>'+INT+')'+SP
    pAtom_occ     = '(?P<Atom_occ>'+DOUBLE+')'+SP
    pAtom_charge  = '(?P<Atom_charge>'+DOUBLE+')'+SP
    pAtom_C60     = '(?P<Atom_C60>'+NINT+')'+SP
    pAtom_name    = '(?P<Atom_name>'+WORD+')'+SP

    if format=="pdb":
        Atom_Record = pAtom_keyword + pAtom_id + pElement_name + pAtom_id1 + pX_val + pY_val + pZ_val + pAtom_charge 
    elif format=="pdb_1":
        Atom_Record = pAtom_keyword + pAtom_id + pElement_name + pAtom_id1 + pX_val + pY_val + pZ_val + pAtom_chain
    elif format=="true_pdb":
        Atom_Record = pAtom_keyword + pAtom_id + pElement_name + pAtom_mol + pAtom_chain + pAtom_id1 + pX_val + pY_val + pZ_val + pAtom_occ + pAtom_charge 
    elif format=="true_pdb2":
        Atom_Record = pAtom_keyword1 + pAtom_id + pAtom_element + pAtom_mol + pAtom_chain + pAtom_id1 + pX_val + pY_val + pZ_val + pAtom_occ + pAtom_charge 
    elif format=="xyz":
        Atom_Record = pElement_name + pX_val + pY_val + pZ_val
    elif format=="iqmol_pdb":
        Atom_Record = pAtom_keyword + pAtom_id + pElement_name + pAtom_chain + pAtom_id1 + pX_val + pY_val + pZ_val + pAtom_occ + pAtom_charge + pAtom_name 




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

            if format=="true_pdb2":
                p1 = '(?P<p1>'+WORD+')'
                s = a[m1.start('Atom_element'):m1.end('Atom_element')]
                m2 = re.search(p1, s)
                if m2!=None:
                    atom_dict["Atom_element"] = s[m2.start('p1'):m2.end('p1')]

            else:
                atom_dict["Atom_element"] = a[m1.start('Element_name'):m1.end('Element_name')]

            atom_dict["Atom_cm_x"] = float(a[m1.start('X_val'):m1.end('X_val')]) * units.Angst
            atom_dict["Atom_cm_y"] = float(a[m1.start('Y_val'):m1.end('Y_val')]) * units.Angst
            atom_dict["Atom_cm_z"] = float(a[m1.start('Z_val'):m1.end('Z_val')]) * units.Angst

            if format=="pdb" or format=="true_pdb":
                atom_dict["Atom_charge"] = float(a[m1.start('Atom_charge'):m1.end('Atom_charge')])
#            atm.Atom_ff_int_type     = int(float(a[m1.start('Atom_type'):m1.end('Atom_type')]))
#            atm.Atom_is_C60_CT       = int(a[m1.start('Atom_C60'):m1.end('Atom_C60')])

            if verbosity>0:
                print("CREATE_ATOM ")
            syst.CREATE_ATOM( Atom(univ,atom_dict)  ) 
    
    #--------- Create bonds ------------------
    print(len(B))
    for b in B:
        m2 = re.search(Bond_Record,b)
        lst = []
        if m2!=None:
            lst  = compINT.findall(m2.group(0))
#            if (int(float(lst[0][0])) == 10):
#               exit(0)
        i = 1   
        while i < len(lst):
            if verbosity>0:
                print('LINK_ATOMS ',lst[0][0], ' and ', lst[i][0],'...')
            syst.LINK_ATOMS(int(float(lst[0][0])),int(float(lst[i][0])))
            i = i + 1

    #----------- Group atoms --------------------
    j = 1
    for fr in C:
        m3 = re.search(Fragment_Record,fr)        
        lst = []
        if m3!=None:
            print(m3.group())
            lst = compINT.findall(m3.group())

            gr_atoms = []
            i = 0
            while i<len(lst):
                gr_atoms.append(int(float(lst[i][0])))
                i = i + 1
#            if len(gr_atoms)>0:
            if verbosity>0:
                print("GROUP_ATOMS",gr_atoms," to form fragment with Group_id = ",j)
            syst.GROUP_ATOMS(gr_atoms,j)
#            syst.show_fragments()
            j = j + 1

    #----------------------------------------------------------------------

