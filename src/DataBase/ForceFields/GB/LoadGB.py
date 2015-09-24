import lcccsObjects
import re



def Load_GB(force_field):

    ff_file = "./DataBase/ForceFields/GB/gb.dat" # In current directory

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
    pFragment_ff_int_type   = '(?P<Fragment_ff_int_type>'+DOUBLE+')'+SP
    pFragment_ff_type       = '(?P<Fragment_ff_type>'+PHRASE+')'+SP
    pFragment_di            = '(?P<Fragment_di>'+DOUBLE+')'+SP
    pFragment_li            = '(?P<Fragment_li>'+DOUBLE+')'+SP
    pFragment_e0            = '(?P<Fragment_e0>'+DOUBLE+')'+SP
    pFragment_rat           = '(?P<Fragment_rat>'+DOUBLE+')'+SP
    pFragment_dw            = '(?P<Fragment_dw>'+DOUBLE+')'+SP
    pFragment_mu            = '(?P<Fragment_mu>'+DOUBLE+')'+SP
    pFragment_nu            = '(?P<Fragment_nu>'+DOUBLE+')'+SP

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


    GB_Fragment_Type_Record = pFragment_ff_int_type + \
                              pFragment_ff_type + \
                              pFragment_di + \
                              pFragment_li + \
                              pFragment_e0 + \
                              pFragment_rat + \
                              pFragment_dw + \
                              pFragment_mu + \
                              pFragment_nu 


    f   = open(ff_file,'r')
    A = f.readlines()
    f.close()


    class atomrecord:
        pass

    ff_par = {}

    for a in A:
        m1 = re.search(GB_Fragment_Type_Record,a)
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
            ff.Fragment_ff_int_type = int(float(a[m1.start('Fragment_ff_int_type'):m1.end('Fragment_ff_int_type')]))
            vFragment_ff_type       = (a[m1.start('Fragment_ff_type'):m1.end('Fragment_ff_type')])
            ff.Fragment_di          = float(a[m1.start('Fragment_di'):m1.end('Fragment_di')])
            ff.Fragment_li          = float(a[m1.start('Fragment_li'):m1.end('Fragment_li')])
            ff.Fragment_e0          = float(a[m1.start('Fragment_e0'):m1.end('Fragment_e0')])
            ff.Fragment_rat         = float(a[m1.start('Fragment_rat'):m1.end('Fragment_rat')])
            ff.Fragment_dw          = float(a[m1.start('Fragment_dw'):m1.end('Fragment_dw')])
            ff.Fragment_mu          = float(a[m1.start('Fragment_mu'):m1.end('Fragment_mu')])
            ff.Fragment_nu          = float(a[m1.start('Fragment_nu'):m1.end('Fragment_nu')])

            # Remove "" from atom type
            sz = len(vAtom_ff_type)     
            ff.Fragment_ff_type = vAtom_ff_type[1:sz-1]

            frag_record = lcccsObjects.Fragment_Record()
            frag_record.set(ff)
            res = force_field.Add_Fragment_Record(frag_record)
            
            print "load type", ff.Fragment_ff_type

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

