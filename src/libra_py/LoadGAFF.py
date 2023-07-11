#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file LoadGAFF.py 
# This module implements functions for loading data into a ForceField object

import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



def Load_GAFF(force_field,  ff_file1 = "data/force_fields/gaff/gaff_types.dat",
                            ff_file2 = "data/force_fields/gaff/gaff_bonds.dat",
                            ff_file3 = "data/force_fields/gaff/gaff_angles.dat",
                            ff_file4 = "data/force_fields/gaff/gaff_dihedrals.dat" ):


#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'
    #DOUBLE = '([-+]?(\d*\.\d*)([eE][-+]?\d+)?)'
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '"((\w|\W)+)"'
    CONNECTOR = '\s*[-]\s*'

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_type_keyword  = '(?P<Atom_type_keyword>'+'TYPE'+')'+SP
    pAtom_ff_int_type   = '(?P<Atom_ff_int_type>'+DOUBLE+')'+SP
    pAtom_ff_type       = '(?P<Atom_ff_type>'+ID+')'+SP
    pAtom1_ff_type      = '(?P<Atom1_ff_type>'+ID+')'
    pAtom2_ff_type      = '(?P<Atom2_ff_type>'+ID+')'
    pAtom3_ff_type      = '(?P<Atom3_ff_type>'+ID+')'
    pAtom4_ff_type      = '(?P<Atom4_ff_type>'+ID+')'

    pBond_r_eq          = '(?P<Bond_r_eq>'+DOUBLE+')'+SP
    pBond_k_bond        = '(?P<Bond_k_bond>'+DOUBLE+')'+SP
    pAngle_theta_eq     = '(?P<Angle_theta_eq>'+DOUBLE+')'+SP
    pAngle_k_angle      = '(?P<Angle_k_angle>'+DOUBLE+')'+SP
    pDihedral_mult      = '(?P<Dihedral_mult>'+INT+')'+SP
    pDihedral_vphi      = '(?P<Dihedral_vphi>'+DOUBLE+')'+SP
    pDihedral_phase     = '(?P<Dihedral_phase>'+DOUBLE+')'+SP
    pAtom_sigma         = '(?P<Atom_sigma>'+DOUBLE+')'+SP
    pAtom_epsilon       = '(?P<Atom_epsilon>'+DOUBLE+')'+SP
    pAtom_C_star        = '(?P<Atom_C_star>'+DOUBLE+')'+SP
    pAtom_Z_star        = '(?P<Atom_Z_star>'+DOUBLE+')'+SP
    pAtom_scale         = '(?P<Atom_scale>'+DOUBLE+')'+SP
    pAtom_GMP           = '(?P<Atom_GMP>'+DOUBLE+')'+SP

    pFF_name_keyword          = '(?P<FF_name_keyword>'+'FORCE_FIELD'+')'+SP
    pFF_sigma_rule_keyword    = '(?P<FF_sigma_rule_keyword>'+'SIGMA_COMB_RULE'+')'+SP
    pFF_epsil_rule_keyword    = '(?P<FF_epsil_rule_keyword>'+'EPSILON_COMB_RULE'+')'+SP
    pFF_vdw_scale13_keyword   = '(?P<FF_vdw_scale13_keyword>'+'VDW_SCALE13'+')'+SP
    pFF_vdw_scale14_keyword   = '(?P<FF_vdw_scale14_keyword>'+'VDW_SCALE14'+')'+SP
    pFF_elec_scale13_keyword  = '(?P<FF_elec_scale13_keyword>'+'ELEC_SCALE13'+')'+SP
    pFF_elec_scale14_keyword  = '(?P<FF_elec_scale14_keyword>'+'ELEC_SCALE14'+')'+SP

    pFF_name_value            = '(?P<FF_name_value>'+WORD+')'+SP
    pFF_sigma_rule_value      = '(?P<FF_sigma_rule_value>'+WORD+')'+SP
    pFF_epsil_rule_value      = '(?P<FF_epsil_rule_value>'+WORD+')'+SP
    pFF_vdw_scale13_value     = '(?P<FF_vdw_scale13_value>'+DOUBLE+')'+SP
    pFF_vdw_scale14_value     = '(?P<FF_vdw_scale14_value>'+DOUBLE+')'+SP
    pFF_elec_scale13_value    = '(?P<FF_elec_scale13_value>'+DOUBLE+')'+SP
    pFF_elec_scale14_value    = '(?P<FF_elec_scale14_value>'+DOUBLE+')'+SP

    FF_name_record = pFF_name_keyword +  pFF_name_value
    FF_sigma_rule_record = pFF_sigma_rule_keyword +  pFF_sigma_rule_value
    FF_epsil_rule_record = pFF_epsil_rule_keyword +  pFF_epsil_rule_value
    FF_vdw_scale13_record = pFF_vdw_scale13_keyword +  pFF_vdw_scale13_value
    FF_vdw_scale14_record = pFF_vdw_scale14_keyword +  pFF_vdw_scale14_value
    FF_elec_scale13_record = pFF_elec_scale13_keyword +  pFF_elec_scale13_value
    FF_elec_scale14_record = pFF_elec_scale14_keyword +  pFF_elec_scale14_value


    GAFF_Atom_Type_Record = pAtom_type_keyword + pAtom_ff_int_type + pAtom_ff_type + pAtom_sigma + pAtom_epsilon + pAtom_C_star + pAtom_Z_star
    GAFF_Bond_Type_Record = pAtom1_ff_type + CONNECTOR + pAtom2_ff_type + SP + pBond_k_bond + pBond_r_eq
    GAFF_Angle_Type_Record = pAtom1_ff_type + CONNECTOR + pAtom2_ff_type + CONNECTOR + pAtom3_ff_type + SP + pAngle_k_angle + pAngle_theta_eq
    GAFF_Dihedral_Type_Record = pAtom1_ff_type + CONNECTOR + pAtom2_ff_type + CONNECTOR + pAtom3_ff_type + CONNECTOR + pAtom4_ff_type + SP + pDihedral_mult + pDihedral_vphi + pDihedral_phase



#-------------- Atoms ------------------------------------
    print("Loading atom types...\n")
    f   = open(ff_file1,'r')
    A = f.readlines()
    f.close()
    class atomrecord:
        pass
    ff_par = {}

    for a in A:
        m1 = re.search(GAFF_Atom_Type_Record,a)
        m2 = re.search(FF_sigma_rule_record,a)
        m3 = re.search(FF_epsil_rule_record,a)
        m4 = re.search(FF_vdw_scale13_record,a)
        m5 = re.search(FF_vdw_scale14_record,a)
        m6 = re.search(FF_elec_scale13_record,a)
        m7 = re.search(FF_elec_scale14_record,a)
        m8 = re.search(FF_name_record,a)

        if m1!=None:
            ff = atomrecord()
            ff.Atom_ff_int_type = int(float(a[m1.start('Atom_ff_int_type'):m1.end('Atom_ff_int_type')]))
            ff.Atom_ff_type     = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            ff.Atom_sigma       = float(a[m1.start('Atom_sigma'):m1.end('Atom_sigma')])
            ff.Atom_epsilon     = float(a[m1.start('Atom_epsilon'):m1.end('Atom_epsilon')])
            ff.Atom_C_star      = float(a[m1.start('Atom_C_star'):m1.end('Atom_C_star')])
            ff.Atom_Z_star      = float(a[m1.start('Atom_Z_star'):m1.end('Atom_Z_star')])          

            atom_record = Atom_Record()
            atom_record.set(ff)
            res = force_field.Add_Atom_Record(atom_record)
            
            print("load type", ff.Atom_ff_type)

        if m2!=None:
#            ff = atomrecord()
            ff_par["sigma_comb_rule"] = a[m2.start('FF_sigma_rule_value'):m2.end('FF_sigma_rule_value')]
#            force_field.set(ff)
        if m3!=None:
#            ff = atomrecord()
            ff_par["epsilon_comb_rule"] = a[m3.start('FF_epsil_rule_value'):m3.end('FF_epsil_rule_value')]
#            force_field.set(ff)
        if m4!=None:
#            ff = atomrecord()
            ff_par["vdw_scale13"] = float(a[m4.start('FF_vdw_scale13_value'):m4.end('FF_vdw_scale13_value')])
#            force_field.set(ff)
        if m5!=None:
#            ff = atomrecord()
            ff_par["vdw_scale14"] = float(a[m5.start('FF_vdw_scale14_value'):m5.end('FF_vdw_scale14_value')])
#            force_field.set(ff)
        if m6!=None:
#            ff = atomrecord()
            ff_par["elec_scale13"] = float(a[m6.start('FF_elec_scale13_value'):m6.end('FF_elec_scale13_value')])
#            force_field.set(ff)
        if m7!=None:
#            ff = atomrecord()
            ff_par["elec_scale14"] = float(a[m7.start('FF_elec_scale14_value'):m7.end('FF_elec_scale14_value')])
#            force_field.set(ff)
        if m8!=None:
#            ff = atomrecord()
            ff_par["ForceField_Name"] = a[m8.start('FF_name_value'):m8.end('FF_name_value')]
#            force_field.set(ff)

    force_field.set(ff_par)




#-------------- Bonds ------------------------------------
    print("Loading bonds...\n")
    f   = open(ff_file2,'r')
    A = f.readlines()
    f.close()
    class bondrecord:
        pass

    for a in A:
        m1 = re.search(GAFF_Bond_Type_Record,a)
        if m1!=None:
            ff = bondrecord()            
            ff.Atom1_ff_type    = a[m1.start('Atom1_ff_type'):m1.end('Atom1_ff_type')]
            ff.Atom2_ff_type    = a[m1.start('Atom2_ff_type'):m1.end('Atom2_ff_type')]
            ff.Bond_k_bond      = float(a[m1.start('Bond_k_bond'):m1.end('Bond_k_bond')])
            ff.Bond_r_eq        = float(a[m1.start('Bond_r_eq'):m1.end('Bond_r_eq')])


            bond_record = Bond_Record()
            bond_record.set(ff)
            #bond_record.show_info()
            res = force_field.Add_Bond_Record(bond_record)            


#----------------- Angles -------------------------------------
    print("Loading angles...\n")
    f   = open(ff_file3,'r')
    A = f.readlines()
    f.close()
    class anglerecord:
        pass

    for a in A:
        m1 = re.search(GAFF_Angle_Type_Record,a)
        if m1!=None:
            ff = anglerecord()
            ff.Atom1_ff_type    = a[m1.start('Atom1_ff_type'):m1.end('Atom1_ff_type')]
            ff.Atom2_ff_type    = a[m1.start('Atom2_ff_type'):m1.end('Atom2_ff_type')]
            ff.Atom3_ff_type    = a[m1.start('Atom3_ff_type'):m1.end('Atom3_ff_type')]
            ff.Angle_k_angle    = float(a[m1.start('Angle_k_angle'):m1.end('Angle_k_angle')])
            ff.Angle_theta_eq   = float(a[m1.start('Angle_theta_eq'):m1.end('Angle_theta_eq')])

            angle_record = Angle_Record()
            angle_record.set(ff)
            #angle_record.show_info()
            res = force_field.Add_Angle_Record(angle_record)


#----------------- Dihedrals -------------------------------------
    print("Loading dihedral...\n")
    f   = open(ff_file4,'r')
    A = f.readlines()
    f.close()
    class dihedralrecord:
        pass

    for a in A:
        m1 = re.search(GAFF_Dihedral_Type_Record,a)
        if m1!=None:
            ff = dihedralrecord()
            ff.Atom1_ff_type    = a[m1.start('Atom1_ff_type'):m1.end('Atom1_ff_type')]
            ff.Atom2_ff_type    = a[m1.start('Atom2_ff_type'):m1.end('Atom2_ff_type')]
            ff.Atom3_ff_type    = a[m1.start('Atom3_ff_type'):m1.end('Atom3_ff_type')]
            ff.Atom4_ff_type    = a[m1.start('Atom4_ff_type'):m1.end('Atom4_ff_type')]
            ff.Dihedral_mult    = int(float(a[m1.start('Dihedral_mult'):m1.end('Dihedral_mult')]))
            ff.Dihedral_vphi    = float(a[m1.start('Dihedral_vphi'):m1.end('Dihedral_vphi')])
            ff.Dihedral_phase   = float(a[m1.start('Dihedral_phase'):m1.end('Dihedral_phase')])

            dihedral_record = Dihedral_Record()
            dihedral_record.set(ff)
            #dihedral_record.show_info()
            res = force_field.Add_Dihedral_Record(dihedral_record)






