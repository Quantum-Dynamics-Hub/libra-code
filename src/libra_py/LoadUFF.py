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
## \file LoadUFF.py 
# This module implements functions for loading data into a ForceField object

import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



def Load_UFF(force_field, ff_file="data/force_fields/uff/uff.dat"):
##
# This function loads data into the force field object provided, assuming a specific format of the input file 
# In this case we assume the input file is formatted to provide some data to set UFF force field
#
# \param[in, out] force_field This is the object which we want to setup
# \param[in] ff_file The name of the file containing all necessary parameters. Default = "uff.dat"
#

#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'
    #DOUBLE = '([-+]?(\d*\.\d*)([eE][-+]?\d+)?)'
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '"((\w|\W)+)"'

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_ff_int_type   = '(?P<Atom_ff_int_type>'+DOUBLE+')'+SP
    pAtom_ff_type       = '(?P<Atom_ff_type>'+PHRASE+')'+SP
    pAtom_radius        = '(?P<Atom_radius>'+DOUBLE+')'+SP
    pAtom_theta         = '(?P<Atom_theta>'+DOUBLE+')'+SP
    pAtom_sigma         = '(?P<Atom_sigma>'+DOUBLE+')'+SP
    pAtom_epsilon       = '(?P<Atom_epsilon>'+DOUBLE+')'+SP
    pAtom_scale         = '(?P<Atom_scale>'+DOUBLE+')'+SP
    pAtom_Z_star        = '(?P<Atom_Z_star>'+DOUBLE+')'+SP
    pAtom_GMP           = '(?P<Atom_GMP>'+DOUBLE+')'+SP

    pFF_name_keyword          = '(?P<FF_name_keyword>'+'FORCE_FIELD'+')'+SP
    pFF_sigma_rule_keyword    = '(?P<FF_sigma_rule_keyword>'+'SIGMA_COMB_RULE'+')'+SP
    pFF_epsil_rule_keyword    = '(?P<FF_epsil_rule_keyword>'+'EPSILON_COMB_RULE'+')'+SP
    pFF_vdw_scale12_keyword   = '(?P<FF_vdw_scale12_keyword>'+'VDW_SCALE12'+')'+SP
    pFF_vdw_scale13_keyword   = '(?P<FF_vdw_scale13_keyword>'+'VDW_SCALE13'+')'+SP
    pFF_vdw_scale14_keyword   = '(?P<FF_vdw_scale14_keyword>'+'VDW_SCALE14'+')'+SP
    pFF_elec_scale12_keyword  = '(?P<FF_elec_scale12_keyword>'+'ELEC_SCALE12'+')'+SP
    pFF_elec_scale13_keyword  = '(?P<FF_elec_scale13_keyword>'+'ELEC_SCALE13'+')'+SP
    pFF_elec_scale14_keyword  = '(?P<FF_elec_scale14_keyword>'+'ELEC_SCALE14'+')'+SP

    pFF_name_value            = '(?P<FF_name_value>'+WORD+')'+SP
    pFF_sigma_rule_value      = '(?P<FF_sigma_rule_value>'+WORD+')'+SP
    pFF_epsil_rule_value      = '(?P<FF_epsil_rule_value>'+WORD+')'+SP
    pFF_vdw_scale12_value     = '(?P<FF_vdw_scale12_value>'+DOUBLE+')'+SP
    pFF_vdw_scale13_value     = '(?P<FF_vdw_scale13_value>'+DOUBLE+')'+SP
    pFF_vdw_scale14_value     = '(?P<FF_vdw_scale14_value>'+DOUBLE+')'+SP
    pFF_elec_scale12_value    = '(?P<FF_elec_scale12_value>'+DOUBLE+')'+SP
    pFF_elec_scale13_value    = '(?P<FF_elec_scale13_value>'+DOUBLE+')'+SP
    pFF_elec_scale14_value    = '(?P<FF_elec_scale14_value>'+DOUBLE+')'+SP

    FF_name_record = pFF_name_keyword +  pFF_name_value
    FF_sigma_rule_record = pFF_sigma_rule_keyword +  pFF_sigma_rule_value
    FF_epsil_rule_record = pFF_epsil_rule_keyword +  pFF_epsil_rule_value
    FF_vdw_scale12_record = pFF_vdw_scale12_keyword +  pFF_vdw_scale12_value
    FF_vdw_scale13_record = pFF_vdw_scale13_keyword +  pFF_vdw_scale13_value
    FF_vdw_scale14_record = pFF_vdw_scale14_keyword +  pFF_vdw_scale14_value
    FF_elec_scale12_record = pFF_elec_scale12_keyword +  pFF_elec_scale12_value
    FF_elec_scale13_record = pFF_elec_scale13_keyword +  pFF_elec_scale13_value
    FF_elec_scale14_record = pFF_elec_scale14_keyword +  pFF_elec_scale14_value


    UFF_Atom_Type_Record = pAtom_ff_int_type + \
                           pAtom_ff_type + \
                           pAtom_radius + \
                           pAtom_theta + \
                           pAtom_sigma + \
                           pAtom_epsilon + \
                           pAtom_scale + \
                           pAtom_Z_star + \
                           pAtom_GMP


    f   = open(ff_file,'r')
    A = f.readlines()
    f.close()


    class atomrecord:
        pass

    ff_par = {}

    for a in A:
        m1 = re.search(UFF_Atom_Type_Record,a)
        m2 = re.search(FF_sigma_rule_record,a)
        m3 = re.search(FF_epsil_rule_record,a)
        m4 = re.search(FF_vdw_scale13_record,a)
        m5 = re.search(FF_vdw_scale14_record,a)
        m6 = re.search(FF_elec_scale13_record,a)
        m7 = re.search(FF_elec_scale14_record,a)
        m8 = re.search(FF_name_record,a)
        m9 = re.search(FF_vdw_scale12_record,a)
        m10= re.search(FF_elec_scale12_record,a)



        if m1!=None:
            ff = atomrecord()
            ff.Atom_ff_int_type = int(float(a[m1.start('Atom_ff_int_type'):m1.end('Atom_ff_int_type')]))
            vAtom_ff_type       = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            ff.Atom_radius      = float(a[m1.start('Atom_radius'):m1.end('Atom_radius')])
            ff.Atom_theta       = float(a[m1.start('Atom_theta'):m1.end('Atom_theta')])
            ff.Atom_sigma       = float(a[m1.start('Atom_sigma'):m1.end('Atom_sigma')])
            ff.Atom_epsilon     = float(a[m1.start('Atom_epsilon'):m1.end('Atom_epsilon')])
            ff.Atom_scale       = float(a[m1.start('Atom_scale'):m1.end('Atom_scale')])
            ff.Atom_Z_star      = float(a[m1.start('Atom_Z_star'):m1.end('Atom_Z_star')])
            ff.Atom_GMP         = float(a[m1.start('Atom_GMP'):m1.end('Atom_GMP')])

            # Remove "" from atom type
            sz = len(vAtom_ff_type)     
            ff.Atom_ff_type = vAtom_ff_type[1:sz-1]

            atom_record = Atom_Record()
            atom_record.set(ff)
            res = force_field.Add_Atom_Record(atom_record)
            
            #print "load type", ff.Atom_ff_type

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
        if m9!=None:
            ff_par["vdw_scale12"] = float(a[m9.start('FF_vdw_scale12_value'):m9.end('FF_vdw_scale12_value')])
        if m10!=None:
            ff_par["elec_scale12"] = float(a[m10.start('FF_elec_scale12_value'):m10.end('FF_elec_scale12_value')])


    force_field.set(ff_par)

