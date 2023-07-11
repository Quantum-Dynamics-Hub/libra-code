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
## \file LoadMMFF94.py 
# This module implements functions for loading data into a ForceField object

import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *




def Load_MMFF94(force_field, 
                ff_file11  = "data/force_fields/mmff94/mmff94_types1.dat",
                ff_file12  = "data/force_fields/mmff94/mmff94_types2.dat",
                ff_file21  = "data/force_fields/mmff94/mmff94_bonds1.dat",
                ff_file22  = "data/force_fields/mmff94/mmff94_bonds2.dat",
                ff_file23  = "data/force_fields/mmff94/mmff94_bonds3.dat",
                ff_file31  = "data/force_fields/mmff94/mmff94_angles1.dat",
                ff_file32  = "data/force_fields/mmff94/mmff94_angles2.dat",
                ff_file33  = "data/force_fields/mmff94/mmff94_angles3.dat",
                ff_file4  = "data/force_fields/mmff94/mmff94_torsions.dat",
                ff_file5  = "data/force_fields/mmff94/mmff94_oop.dat" ):


#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'
    #DOUBLE = '([-+]?(\d*\.\d*)([eE][-+]?\d+)?)'
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    WORD1  = '(([a-zA-Z]|[-])+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '"((\w|\W)+)"'
    CONNECTOR = '\s*[-]\s*'

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_type_keyword  = '(?P<Atom_type_keyword>'+'TYPE'+')'+SP
    pAtom_ff_int_type   = '(?P<Atom_ff_int_type>'+DOUBLE+')'+SP
    pAtom_ff_type       = '(?P<Atom_ff_type>'+PHRASE+')'+SP
    pAtom_ff_type_H     = '(?P<Atom_ff_type_H>'+PHRASE+')'+SP
    pAtom_atomic_number = '(?P<Atom_atomic_number>'+DOUBLE+')'+SP
    pAtom1_ff_type      = '(?P<Atom1_ff_type>'+ID+')'+SP
    pAtom2_ff_type      = '(?P<Atom2_ff_type>'+ID+')'+SP
    pAtom3_ff_type      = '(?P<Atom3_ff_type>'+ID+')'+SP
    pAtom4_ff_type      = '(?P<Atom4_ff_type>'+ID+')'+SP
    pAtom1_ff_int_type  = '(?P<Atom1_ff_int_type>'+DOUBLE+')'+SP
    pAtom2_ff_int_type  = '(?P<Atom2_ff_int_type>'+DOUBLE+')'+SP
    pAtom3_ff_int_type  = '(?P<Atom3_ff_int_type>'+DOUBLE+')'+SP
    pAtom4_ff_int_type  = '(?P<Atom4_ff_int_type>'+DOUBLE+')'+SP


    pX_type_index       = '(?P<X_type_index>'+DOUBLE+')'+SP
    pBond_r_eq          = '(?P<Bond_r_eq>'+DOUBLE+')'+SP
    pBond_k_bond        = '(?P<Bond_k_bond>'+DOUBLE+')'+SP
    pBond_bci           = '(?P<Bond_bci>'+DOUBLE+')'+SP
    pAngle_theta_eq     = '(?P<Angle_theta_eq>'+DOUBLE+')'+SP
    pAngle_k_angle      = '(?P<Angle_k_angle>'+DOUBLE+')'+SP
    pAngle_k_ijk        = '(?P<Angle_k_ijk>'+DOUBLE+')'+SP
    pAngle_k_kji        = '(?P<Angle_k_kji>'+DOUBLE+')'+SP
    pDihedral_mult      = '(?P<Dihedral_mult>'+INT+')'+SP
    pDihedral_vphi      = '(?P<Dihedral_vphi>'+DOUBLE+')'+SP
    pDihedral_phase     = '(?P<Dihedral_phase>'+DOUBLE+')'+SP
    pTors_V1            = '(?P<Tors_V1>'+DOUBLE+')'+SP
    pTors_V2            = '(?P<Tors_V2>'+DOUBLE+')'+SP
    pTors_V3            = '(?P<Tors_V3>'+DOUBLE+')'+SP




    pAtom_crd           = '(?P<Atom_crd>'+DOUBLE+')'+SP
    pAtom_val           = '(?P<Atom_val>'+DOUBLE+')'+SP
    pAtom_pilp          = '(?P<Atom_pilp>'+DOUBLE+')'+SP
    pAtom_mltb          = '(?P<Atom_mltb>'+DOUBLE+')'+SP
    pAtom_arom          = '(?P<Atom_arom>'+DOUBLE+')'+SP
    pAtom_lin           = '(?P<Atom_lin>'+DOUBLE+')'+SP
    pAtom_sbmb          = '(?P<Atom_sbmb>'+DOUBLE+')'+SP
    pAtom_alpha         = '(?P<Atom_alpha>'+DOUBLE+')'+SP
    pAtom_N_eff         = '(?P<Atom_N_eff>'+DOUBLE+')'+SP
    pAtom_A_scale       = '(?P<Atom_A_scale>'+DOUBLE+')'+SP
    pAtom_G_scale       = '(?P<Atom_G_scale>'+DOUBLE+')'+SP
    pAtom_DAN           = '(?P<Atom_DAN>'+WORD+')'+SP
    pAtom_pbci          = '(?P<Atom_pbci>'+DOUBLE+')'+SP
    pAtom_fcadj         = '(?P<Atom_fcadj>'+DOUBLE+')'+SP

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

    pFF_name_value            = '(?P<FF_name_value>'+ID+')'+SP
    pFF_sigma_rule_value      = '(?P<FF_sigma_rule_value>'+WORD1+')'+SP
    pFF_epsil_rule_value      = '(?P<FF_epsil_rule_value>'+WORD1+')'+SP
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


    MMFF94_Atom_Type_Record1 = pAtom_type_keyword + \
                               pAtom_ff_type + \
                               pAtom_ff_type_H + \
                               pAtom_ff_int_type + \
                               pAtom_atomic_number + \
                               pAtom_crd + \
                               pAtom_val + \
                               pAtom_pilp + \
                               pAtom_mltb + \
                               pAtom_arom + \
                               pAtom_lin + \
                               pAtom_sbmb + \
                               pAtom_alpha + \
                               pAtom_N_eff + \
                               pAtom_A_scale + \
                               pAtom_G_scale + \
                               pAtom_DAN + \
                               pAtom_pbci + \
                               pAtom_fcadj

    MMFF94_Atom_Type_Record2 = pAtom_type_keyword + \
                               pAtom_ff_type + \
                               pAtom_ff_int_type + \
                               pAtom1_ff_int_type + \
                               pAtom2_ff_int_type + \
                               pAtom3_ff_int_type + \
                               pAtom4_ff_int_type 
                          
    MMFF94_Bond_Type_Record1 = pX_type_index + pAtom1_ff_int_type + pAtom2_ff_int_type + pBond_k_bond + pBond_r_eq
    MMFF94_Bond_Type_Record2 = pX_type_index + pAtom1_ff_int_type + pAtom2_ff_int_type + pBond_bci
    MMFF94_Bond_Type_Record3 = pAtom1_ff_int_type + pAtom2_ff_int_type + pBond_r_eq + pBond_k_bond

    MMFF94_Angle_Type_Record1 = pX_type_index + pAtom1_ff_int_type + pAtom2_ff_int_type + pAtom3_ff_int_type + pAngle_k_angle + pAngle_theta_eq
    MMFF94_Angle_Type_Record2 = pX_type_index + pAtom1_ff_int_type + pAtom2_ff_int_type + pAtom3_ff_int_type + pAngle_k_ijk + pAngle_k_kji


    MMFF94_Dihedral_Type_Record = pX_type_index + pAtom1_ff_int_type + pAtom2_ff_int_type + pAtom3_ff_int_type + pAtom4_ff_int_type + pTors_V1 + pTors_V2 + pTors_V3

    MMFF94_OOP_Type_Record = pAtom1_ff_int_type + pAtom2_ff_int_type + pAtom3_ff_int_type + pAtom4_ff_int_type + pDihedral_vphi



#-------------- Atoms ------------------------------------
    print("Loading atom types...\n")
    f   = open(ff_file11,'r')
    A1 = f.readlines()
    f.close()
    f   = open(ff_file12,'r')
    A2 = f.readlines()
    f.close()

    class atomrecord:
        pass

    ff_par = {}

    for a in A1:
        m1 = re.search(MMFF94_Atom_Type_Record1,a)
        m2 = re.search(FF_sigma_rule_record,a)
        m3 = re.search(FF_epsil_rule_record,a)
        m4 = re.search(FF_vdw_scale13_record,a)
        m5 = re.search(FF_vdw_scale14_record,a)
        m6 = re.search(FF_elec_scale13_record,a)
        m7 = re.search(FF_elec_scale14_record,a)
        m8 = re.search(FF_name_record,a)

        if m1!=None:
            ff = atomrecord()
            vAtom_ff_type       = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            # Remove "" from atom type
            sz = len(vAtom_ff_type)
            ff.Atom_ff_type = vAtom_ff_type[1:sz-1]
            
            vAtom_ff_type_H     = (a[m1.start('Atom_ff_type_H'):m1.end('Atom_ff_type_H')])
            # Remove "" from atom type H
            sz = len(vAtom_ff_type_H)
            ff.Atom_ff_type_H = vAtom_ff_type_H[1:sz-1]
            
            ff.Atom_ff_int_type = int(float(a[m1.start('Atom_ff_int_type'):m1.end('Atom_ff_int_type')]))
            ff.Atom_atomic_number = int(float(a[m1.start('Atom_atomic_number'):m1.end('Atom_atomic_number')]))
            ff.Atom_crd         = int(float(a[m1.start('Atom_crd'):m1.end('Atom_crd')]))
            ff.Atom_val         = int(float(a[m1.start('Atom_val'):m1.end('Atom_val')]))
            ff.Atom_pilp        = int(float(a[m1.start('Atom_pilp'):m1.end('Atom_pilp')]))
            ff.Atom_mltb        = int(float(a[m1.start('Atom_mltb'):m1.end('Atom_mltb')]))
            ff.Atom_arom        = int(float(a[m1.start('Atom_arom'):m1.end('Atom_arom')]))
            ff.Atom_lin         = int(float(a[m1.start('Atom_lin'):m1.end('Atom_lin')]))
            ff.Atom_sbmb        = int(float(a[m1.start('Atom_sbmb'):m1.end('Atom_sbmb')]))
            ff.Atom_alpha       = float(a[m1.start('Atom_alpha'):m1.end('Atom_alpha')])
            ff.Atom_N_eff       = float(a[m1.start('Atom_N_eff'):m1.end('Atom_N_eff')])
            ff.Atom_A_scale     = float(a[m1.start('Atom_A_scale'):m1.end('Atom_A_scale')])
            ff.Atom_G_scale     = float(a[m1.start('Atom_G_scale'):m1.end('Atom_G_scale')])
            ff.Atom_DAN         = a[m1.start('Atom_DAN'):m1.end('Atom_DAN')]
            ff.Atom_pbci        = float(a[m1.start('Atom_pbci'):m1.end('Atom_pbci')])
            ff.Atom_fcadj       = float(a[m1.start('Atom_fcadj'):m1.end('Atom_fcadj')])

            atom_record = Atom_Record()
            atom_record.set(ff)
            res = force_field.Add_Atom_Record(atom_record)                       
            print("load type", ff.Atom_ff_type)

        if m2!=None:
            ff_par["sigma_comb_rule"] = a[m2.start('FF_sigma_rule_value'):m2.end('FF_sigma_rule_value')]
        if m3!=None:
            ff_par["epsilon_comb_rule"] = a[m3.start('FF_epsil_rule_value'):m3.end('FF_epsil_rule_value')]
        if m4!=None:
            ff_par["vdw_scale13"] = float(a[m4.start('FF_vdw_scale13_value'):m4.end('FF_vdw_scale13_value')])
        if m5!=None:
            ff_par["vdw_scale14"] = float(a[m5.start('FF_vdw_scale14_value'):m5.end('FF_vdw_scale14_value')])
        if m6!=None:
            ff_par["elec_scale13"] = float(a[m6.start('FF_elec_scale13_value'):m6.end('FF_elec_scale13_value')])
        if m7!=None:
            ff_par["elec_scale14"] = float(a[m7.start('FF_elec_scale14_value'):m7.end('FF_elec_scale14_value')])
        if m8!=None:
            ff_par["ForceField_Name"] = a[m8.start('FF_name_value'):m8.end('FF_name_value')]

    force_field.set(ff_par)

    for a in A2:
        m1 = re.search(MMFF94_Atom_Type_Record2,a)

        if m1!=None:
            ff = atomrecord()
            vAtom_ff_type       = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            # Remove "" from atom type
            sz = len(vAtom_ff_type)
            ff.Atom_ff_type = vAtom_ff_type[1:sz-1]

            ff.Atom_ff_int_type = int(float(a[m1.start('Atom_ff_int_type'):m1.end('Atom_ff_int_type')]))
            ff.Atom_ff_eq_int_type2 = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom_ff_eq_int_type3 = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Atom_ff_eq_int_type4 = int(float(a[m1.start('Atom3_ff_int_type'):m1.end('Atom3_ff_int_type')]))
            ff.Atom_ff_eq_int_type5 = int(float(a[m1.start('Atom4_ff_int_type'):m1.end('Atom4_ff_int_type')]))

            atom_record = Atom_Record()
            atom_record.set(ff)
            res = force_field.Add_Atom_Record(atom_record)


 

#    force_field.show_atom_records()
#    force_field.show_info()


#-------------- Bonds ------------------------------------
    print("Loading bonds...\n")
    f   = open(ff_file21,'r')
    A1 = f.readlines()
    f.close()
    f   = open(ff_file22,'r')
    A2 = f.readlines()
    f.close()
    f   = open(ff_file23,'r')
    A3 = f.readlines()
    f.close()

    class bondrecord:
        pass

    for a in A1:
        m1 = re.search(MMFF94_Bond_Type_Record1,a)
        if m1!=None:
            ff = bondrecord()            
            ff.Atom1_ff_int_type = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Bond_type_index  = int(float(a[m1.start('X_type_index'):m1.end('X_type_index')])) 
            ff.Bond_k_bond      = float(a[m1.start('Bond_k_bond'):m1.end('Bond_k_bond')])
            ff.Bond_r_eq        = float(a[m1.start('Bond_r_eq'):m1.end('Bond_r_eq')])

            bond_record = Bond_Record()
            bond_record.set(ff)
            res = force_field.Add_Bond_Record(bond_record)            

    for a in A2:
        m1 = re.search(MMFF94_Bond_Type_Record2,a)
        if m1!=None:
#            print a
            ff = bondrecord()
            ff.Atom1_ff_int_type = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Bond_type_index  = int(float(a[m1.start('X_type_index'):m1.end('X_type_index')]))
            ff.Bond_bci         = float(a[m1.start('Bond_bci'):m1.end('Bond_bci')])         

            bond_record = Bond_Record()
            bond_record.set(ff)
            res = force_field.Add_Bond_Record(bond_record)

    print("Setting parameters for rule-based calculations")
    for a in A3:
        m1 = re.search(MMFF94_Bond_Type_Record3,a)     
        # Note that the meaning of patterns is different from that of variables!!!
        if m1!=None:
            ff = bondrecord()
            ff.Atom1_atomic_number = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_atomic_number = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Bond_r_eq_ref    = float(a[m1.start('Bond_r_eq'):m1.end('Bond_r_eq')])
            ff.Bond_k_bond_ref  = float(a[m1.start('Bond_k_bond'):m1.end('Bond_k_bond')])

            bond_record = Bond_Record()
            bond_record.set(ff)
#            bond_record.show_info()
            res = force_field.Add_Bond_Record(bond_record)
            print(res)



#    print "Printing all bond records in this force field"
#    force_field.show_bond_records();



#----------------- Angles -------------------------------------
    print("Loading angles and stretch-bending...\n")
    f   = open(ff_file31,'r')
    A1 = f.readlines()
    f.close()
    f   = open(ff_file32,'r')
    A2 = f.readlines()
    f.close()
    f   = open(ff_file33,'r')
    A3 = f.readlines()
    f.close()



    class anglerecord:
        pass

    for a in A1:
        m1 = re.search(MMFF94_Angle_Type_Record1,a)
        if m1!=None:
            ff = anglerecord()
            ff.Atom1_ff_int_type  = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type  = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Atom3_ff_int_type  = int(float(a[m1.start('Atom3_ff_int_type'):m1.end('Atom3_ff_int_type')]))
            ff.Angle_type_index   = int(float(a[m1.start('X_type_index'):m1.end('X_type_index')]))
            ff.Angle_k_angle      = float(a[m1.start('Angle_k_angle'):m1.end('Angle_k_angle')])
            ff.Angle_theta_eq     = float(a[m1.start('Angle_theta_eq'):m1.end('Angle_theta_eq')])

            angle_record = Angle_Record()
            angle_record.set(ff) 
            res = force_field.Add_Angle_Record(angle_record,0)

    for a in A2:
        m1 = re.search(MMFF94_Angle_Type_Record2,a)
        if m1!=None:
            ff = anglerecord()
            ff.Atom1_ff_int_type  = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type  = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Atom3_ff_int_type  = int(float(a[m1.start('Atom3_ff_int_type'):m1.end('Atom3_ff_int_type')]))
            ff.Angle_type_index   = int(float(a[m1.start('X_type_index'):m1.end('X_type_index')]))
            ff.Angle_kijk_sb      = float(a[m1.start('Angle_k_ijk'):m1.end('Angle_k_ijk')])
            ff.Angle_kkji_sb      = float(a[m1.start('Angle_k_kji'):m1.end('Angle_k_kji')])

            angle_record = Angle_Record()
            angle_record.set(ff)
            res = force_field.Add_Angle_Record(angle_record,1)


#    print "Printing all angle records in this force field"
#    force_field.show_angle_records();



#----------------- Dihedrals -------------------------------------
    print("Loading dihedral...\n")
    f   = open(ff_file4,'r')
    A = f.readlines()
    f.close()
    class dihedralrecord:
        pass

    for a in A:
        m1 = re.search(MMFF94_Dihedral_Type_Record,a)
        if m1!=None:
            ff = dihedralrecord()
            ff.Atom1_ff_int_type  = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type  = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Atom3_ff_int_type  = int(float(a[m1.start('Atom3_ff_int_type'):m1.end('Atom3_ff_int_type')]))
            ff.Atom4_ff_int_type  = int(float(a[m1.start('Atom4_ff_int_type'):m1.end('Atom4_ff_int_type')]))
            ff.Dihedral_type_index= int(float(a[m1.start('X_type_index'):m1.end('X_type_index')]))
            ff.Dihedral_vphi1     = float(a[m1.start('Tors_V1'):m1.end('Tors_V1')])
            ff.Dihedral_vphi2     = float(a[m1.start('Tors_V2'):m1.end('Tors_V2')])
            ff.Dihedral_vphi3     = float(a[m1.start('Tors_V3'):m1.end('Tors_V3')])
            ff.Dihedral_vphi      = float(a[m1.start('Tors_V1'):m1.end('Tors_V1')])


            dihedral_record = Dihedral_Record()
            dihedral_record.set(ff)
            res = force_field.Add_Dihedral_Record(dihedral_record,0)

#    print "Printing all dihedral(torsion) records in this force field"
#    force_field.show_dihedral_records();


#----------------- Out-of-plane -------------------------------------
    print("Loading oop...\n")
    f   = open(ff_file5,'r')
    A = f.readlines()
    f.close()
    class dihedralrecord:
        pass

    for a in A:
        m1 = re.search(MMFF94_OOP_Type_Record,a)
        if m1!=None:
            ff = dihedralrecord()
            ff.Atom1_ff_int_type  = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type  = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Atom3_ff_int_type  = int(float(a[m1.start('Atom3_ff_int_type'):m1.end('Atom3_ff_int_type')]))
            ff.Atom4_ff_int_type  = int(float(a[m1.start('Atom4_ff_int_type'):m1.end('Atom4_ff_int_type')]))

            ff.Dihedral_vphi      = float(a[m1.start('Dihedral_vphi'):m1.end('Dihedral_vphi')])


            dihedral_record = Dihedral_Record()
            dihedral_record.set(ff)
            res = force_field.Add_Improper_Record(dihedral_record)

#    print "Printing all improper(oop) records in this force field"
#    force_field.show_improper_records();






