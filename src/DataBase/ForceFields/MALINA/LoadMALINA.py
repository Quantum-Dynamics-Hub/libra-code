import lcccsObjects
import re



def Load_MALINA(force_field):

    ff_file1  = "./DataBase/ForceFields/MALINA/malina_types.dat" # In current directory
    ff_file2  = "./DataBase/ForceFields/MALINA/malina_bonds1.dat" # In current directory

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
    pAtom_charge        = '(?P<Atom_charge>'+DOUBLE+')'+SP
    pAtom1_ff_type      = '(?P<Atom1_ff_type>'+ID+')'+SP
    pAtom2_ff_type      = '(?P<Atom2_ff_type>'+ID+')'+SP
    pAtom3_ff_type      = '(?P<Atom3_ff_type>'+ID+')'+SP
    pAtom4_ff_type      = '(?P<Atom4_ff_type>'+ID+')'+SP
    pAtom1_ff_int_type  = '(?P<Atom1_ff_int_type>'+DOUBLE+')'+SP
    pAtom2_ff_int_type  = '(?P<Atom2_ff_int_type>'+DOUBLE+')'+SP
    pAtom3_ff_int_type  = '(?P<Atom3_ff_int_type>'+DOUBLE+')'+SP
    pAtom4_ff_int_type  = '(?P<Atom4_ff_int_type>'+DOUBLE+')'+SP

    pBond_r_eq          = '(?P<Bond_r_eq>'+DOUBLE+')'+SP
    pBond_D_bond        = '(?P<Bond_D_bond>'+DOUBLE+')'+SP
    pBond_alpha         = '(?P<Bond_alpha>'+DOUBLE+')'+SP
    pBond_delta         = '(?P<Bond_delta>'+DOUBLE+')'+SP
    pBond_bci           = '(?P<Bond_bci>'+DOUBLE+')'+SP
    pBond_wij           = '(?P<Bond_wij>'+DOUBLE+')'+SP
    pBond_wij_1         = '(?P<Bond_wij_1>'+DOUBLE+')'+SP
    pBond_wij_2         = '(?P<Bond_wij_2>'+DOUBLE+')'+SP
    pBond_alpij         = '(?P<Bond_alpij>'+DOUBLE+')'+SP
    pBond_alpij_1       = '(?P<Bond_alpij_1>'+DOUBLE+')'+SP
    pBond_alpij_2       = '(?P<Bond_alpij_2>'+DOUBLE+')'+SP
    pBond_rij_1         = '(?P<Bond_rij_1>'+DOUBLE+')'+SP
    pBond_rij_2         = '(?P<Bond_rij_2>'+DOUBLE+')'+SP




    pFF_sigma_rule_keyword    = '(?P<FF_sigma_rule_keyword>'+'SIGMA_COMB_RULE'+')'+SP
    pFF_epsil_rule_keyword    = '(?P<FF_epsil_rule_keyword>'+'EPSILON_COMB_RULE'+')'+SP
    pFF_vdw_scale13_keyword   = '(?P<FF_vdw_scale13_keyword>'+'VDW_SCALE13'+')'+SP
    pFF_vdw_scale14_keyword   = '(?P<FF_vdw_scale14_keyword>'+'VDW_SCALE14'+')'+SP
    pFF_elec_scale13_keyword  = '(?P<FF_elec_scale13_keyword>'+'ELEC_SCALE13'+')'+SP
    pFF_elec_scale14_keyword  = '(?P<FF_elec_scale14_keyword>'+'ELEC_SCALE14'+')'+SP

    pFF_sigma_rule_value      = '(?P<FF_sigma_rule_value>'+WORD1+')'+SP
    pFF_epsil_rule_value      = '(?P<FF_epsil_rule_value>'+WORD1+')'+SP
    pFF_vdw_scale13_value     = '(?P<FF_vdw_scale13_value>'+DOUBLE+')'+SP
    pFF_vdw_scale14_value     = '(?P<FF_vdw_scale14_value>'+DOUBLE+')'+SP
    pFF_elec_scale13_value    = '(?P<FF_elec_scale13_value>'+DOUBLE+')'+SP
    pFF_elec_scale14_value    = '(?P<FF_elec_scale14_value>'+DOUBLE+')'+SP

    FF_sigma_rule_record = pFF_sigma_rule_keyword +  pFF_sigma_rule_value
    FF_epsil_rule_record = pFF_epsil_rule_keyword +  pFF_epsil_rule_value
    FF_vdw_scale13_record = pFF_vdw_scale13_keyword +  pFF_vdw_scale13_value
    FF_vdw_scale14_record = pFF_vdw_scale14_keyword +  pFF_vdw_scale14_value
    FF_elec_scale13_record = pFF_elec_scale13_keyword +  pFF_elec_scale13_value
    FF_elec_scale14_record = pFF_elec_scale14_keyword +  pFF_elec_scale14_value


    MALINA_Atom_Type_Record = pAtom_type_keyword + \
                              pAtom_ff_type + \
                              pAtom_ff_int_type + \
                              pAtom_atomic_number + \
                              pAtom_charge

    MALINA_Bond_Type_Record = pAtom1_ff_int_type + \
                              pAtom2_ff_int_type + \
                              pBond_D_bond + \
                              pBond_r_eq + \
                              pBond_alpha + \
                              pBond_delta + \
                              pBond_bci + \
                              pBond_wij_1 + \
                              pBond_wij_2 + \
                              pBond_alpij_1 + \
                              pBond_alpij_2 + \
                              pBond_rij_1 + \
                              pBond_rij_2 + \
                              pBond_wij + \
                              pBond_alpij 


#-------------- Atoms ------------------------------------
    print "Loading atom types...\n"
    f   = open(ff_file1,'r')
    A = f.readlines()
    f.close()
    class atomrecord:
        pass

    for a in A:
        m1 = re.search(MALINA_Atom_Type_Record,a)
        m2 = re.search(FF_sigma_rule_record,a)
        m3 = re.search(FF_epsil_rule_record,a)
        m4 = re.search(FF_vdw_scale13_record,a)
        m5 = re.search(FF_vdw_scale14_record,a)
        m6 = re.search(FF_elec_scale13_record,a)
        m7 = re.search(FF_elec_scale14_record,a)

        if m1!=None:
            ff = atomrecord()
            vAtom_ff_type       = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            # Remove "" from atom type
            sz = len(vAtom_ff_type)
            ff.Atom_ff_type = vAtom_ff_type[1:sz-1]
            
            
            ff.Atom_ff_int_type = int(float(a[m1.start('Atom_ff_int_type'):m1.end('Atom_ff_int_type')]))
            ff.Atom_atomic_number = int(float(a[m1.start('Atom_atomic_number'):m1.end('Atom_atomic_number')]))
            ff.Atom_Z_star = float(a[m1.start('Atom_charge'):m1.end('Atom_charge')])

           
            atom_record = lcccsObjects.Atom_Record()
            atom_record.set(ff)
#            atom_record.show_info()
            res = force_field.Add_Atom_Record(atom_record)                       
            print "load type", ff.Atom_ff_type
#            print res

        if m2!=None:
            ff = atomrecord()
            ff.ForceField_sigma_comb_rule = a[m2.start('FF_sigma_rule_value'):m2.end('FF_sigma_rule_value')]
            force_field.set(ff)
        if m3!=None:
            ff = atomrecord()
            ff.ForceField_epsilon_comb_rule = a[m3.start('FF_epsil_rule_value'):m3.end('FF_epsil_rule_value')]
            force_field.set(ff)
        if m4!=None:
            ff = atomrecord()
            ff.ForceField_vdw_scale13 = float(a[m4.start('FF_vdw_scale13_value'):m4.end('FF_vdw_scale13_value')])
            force_field.set(ff)
        if m5!=None:
            ff = atomrecord()
            ff.ForceField_vdw_scale14 = float(a[m5.start('FF_vdw_scale14_value'):m5.end('FF_vdw_scale14_value')])
            force_field.set(ff)
        if m6!=None:
            ff = atomrecord()
            ff.ForceField_elec_scale13 = float(a[m6.start('FF_elec_scale13_value'):m6.end('FF_elec_scale13_value')])
            force_field.set(ff)
        if m7!=None:
            ff = atomrecord()
            ff.ForceField_elec_scale14 = float(a[m7.start('FF_elec_scale14_value'):m7.end('FF_elec_scale14_value')])
            force_field.set(ff)

    force_field.show_atom_records()
    force_field.show_info()


#-------------- Bonds ------------------------------------
    print "Loading bonds...\n"
    f   = open(ff_file2,'r')
    A1 = f.readlines()
    f.close()

    class bondrecord:
        pass

    for a in A1:
        m1 = re.search(MALINA_Bond_Type_Record,a)
        if m1!=None:
            ff = bondrecord()            
            ff.Atom1_ff_int_type = int(float(a[m1.start('Atom1_ff_int_type'):m1.end('Atom1_ff_int_type')]))
            ff.Atom2_ff_int_type = int(float(a[m1.start('Atom2_ff_int_type'):m1.end('Atom2_ff_int_type')]))
            ff.Bond_D_bond      = float(a[m1.start('Bond_D_bond'):m1.end('Bond_D_bond')])
            ff.Bond_r_eq        = float(a[m1.start('Bond_r_eq'):m1.end('Bond_r_eq')])
            ff.Bond_alpha       = float(a[m1.start('Bond_alpha'):m1.end('Bond_alpha')])
            ff.Bond_shift_elec  = float(a[m1.start('Bond_delta'):m1.end('Bond_delta')])
            ff.Bond_bci         = float(a[m1.start('Bond_bci'):m1.end('Bond_bci')])
            ff.Bond_wij         = float(a[m1.start('Bond_wij'):m1.end('Bond_wij')])
            ff.Bond_wij_1       = float(a[m1.start('Bond_wij_1'):m1.end('Bond_wij_1')])
            ff.Bond_wij_2       = float(a[m1.start('Bond_wij_2'):m1.end('Bond_wij_2')])
            ff.Bond_alpij       = float(a[m1.start('Bond_alpij'):m1.end('Bond_alpij')])
            ff.Bond_alpij_1     = float(a[m1.start('Bond_alpij_1'):m1.end('Bond_alpij_1')])
            ff.Bond_alpij_2     = float(a[m1.start('Bond_alpij_2'):m1.end('Bond_alpij_2')])
            ff.Bond_rij_1       = float(a[m1.start('Bond_rij_1'):m1.end('Bond_rij_1')])
            ff.Bond_rij_2       = float(a[m1.start('Bond_rij_2'):m1.end('Bond_rij_2')])


            bond_record = lcccsObjects.Bond_Record()
            bond_record.set(ff)
            res = force_field.Add_Bond_Record(bond_record)           

    print "Printing all bond records in this force field"
    force_field.show_bond_records();
 






