import lcccsObjects
import re



def Load_ESFF(force_field):

    ff_file = "./DataBase/ForceFields/ESFF/esff1.dat" # In current directory

#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'
    #DOUBLE = '([-+]?(\d*\.\d*)([eE][-+]?\d+)?)'
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '((\w|[-+=*]|\')+)' 

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pAtom_ff_int_type   = '(?P<Atom_ff_int_type>'+DOUBLE+')'+SP
    pAtom_ff_type       = '(?P<Atom_ff_type>'+PHRASE+')'+SP
    pAtom_dative        = '(?P<Atom_dative>'+DOUBLE+')'+SP
    pAtom_brdr1         = '(?P<Atom_brdr1>'+DOUBLE+')'+SP
    pAtom_brdr2         = '(?P<Atom_brdr2>'+DOUBLE+')'+SP
    pAtom_brdr3         = '(?P<Atom_brdr3>'+DOUBLE+')'+SP

    pAtom_theta         = '(?P<Atom_theta>'+DOUBLE+')'+SP
    pAtom_sigma         = '(?P<Atom_sigma>'+DOUBLE+')'+SP
    pAtom_epsilon       = '(?P<Atom_epsilon>'+DOUBLE+')'+SP
    pAtom_scale         = '(?P<Atom_scale>'+DOUBLE+')'+SP
    pAtom_Z_star        = '(?P<Atom_Z_star>'+DOUBLE+')'+SP
    pAtom_GMP           = '(?P<Atom_GMP>'+DOUBLE+')'+SP

    ESFF_Bond_Ref_Value_Record3 = pAtom_ff_type + \
                                 pAtom_dative + \
                                 pAtom_brdr1 + \
                                 pAtom_brdr2 + \
                                 pAtom_brdr3 

    ESFF_Bond_Ref_Value_Record2 = pAtom_ff_type + \
                                 pAtom_dative + \
                                 pAtom_brdr1 + \
                                 pAtom_brdr2 

    ESFF_Bond_Ref_Value_Record1 = pAtom_ff_type + \
                                 pAtom_dative + \
                                 pAtom_brdr1                          
                           



    f   = open(ff_file,'r')
    A = f.readlines()
    f.close()


    class atomrecord:
        pass

    for a in A:
        m1 = re.search(ESFF_Bond_Ref_Value_Record1,a)
        m2 = re.search(ESFF_Bond_Ref_Value_Record2,a)
        m3 = re.search(ESFF_Bond_Ref_Value_Record3,a)        
        if m1!=None:
            ff = atomrecord()
            
            ff.Atom_ff_type     = (a[m1.start('Atom_ff_type'):m1.end('Atom_ff_type')])
            ff.Atom_dative      = float(a[m1.start('Atom_dative'):m1.end('Atom_dative')])
            ff.Atom_brdr1       = float(a[m1.start('Atom_brdr1'):m1.end('Atom_brdr1')])

            if m2!=None:
                ff.Atom_brdr2       = float(a[m2.start('Atom_brdr2'):m2.end('Atom_brdr2')])
            if m3!=None:
                ff.Atom_brdr3       = float(a[m3.start('Atom_brdr3'):m3.end('Atom_brdr3')])
                               

            atom_record = lcccsObjects.Atom_Record()
            atom_record.set(ff)

            res = force_field.Add_Atom_Record(atom_record)
            
            print "load type", ff.Atom_ff_type
#            atom_record.show_info()     


