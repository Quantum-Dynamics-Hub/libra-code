import lcccsObjects
import re



def Load_SuttonChen(force_field):

    ff_file = "./DataBase/ForceFields/SuttonChen/suttonchen.dat"

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
    pAtom_element      = '(?P<Atom_element>'+WORD+')'+SP
    pAtom_such_n       = '(?P<Atom_such_n>'+INT+')'+SP
    pAtom_such_m       = '(?P<Atom_such_m>'+INT+')'+SP
    pAtom_such_a       = '(?P<Atom_such_a>'+DOUBLE+')'+SP
    pAtom_such_D       = '(?P<Atom_such_D>'+DOUBLE+')'+SP
    pAtom_such_c       = '(?P<Atom_such_c>'+DOUBLE+')'+SP   

    SUCH_Atom_Type_Record = pAtom_element + \
                            pAtom_such_n + \
                            pAtom_such_m + \
                            pAtom_such_D + \
                            pAtom_such_c + \
                            pAtom_such_a


    f   = open(ff_file,'r')
    A = f.readlines()
    f.close()


    class atomrecord:
        pass

    convert = 0.0964853 # convertion from meV to kcal/mol
    for a in A:
        m1 = re.search(SUCH_Atom_Type_Record,a)
        if m1!=None:
            ff = atomrecord()
            ff.Atom_ff_type     = a[m1.start('Atom_element'):m1.end('Atom_element')]   
            ff.Atom_element     = a[m1.start('Atom_element'):m1.end('Atom_element')]
            ff.Atom_such_n      = float(a[m1.start('Atom_such_n'):m1.end('Atom_such_n')])
            ff.Atom_such_m      = float(a[m1.start('Atom_such_m'):m1.end('Atom_such_m')])
            ff.Atom_such_a      = float(a[m1.start('Atom_such_a'):m1.end('Atom_such_a')])
            ff.Atom_such_D      = convert*float(a[m1.start('Atom_such_D'):m1.end('Atom_such_D')])
            ff.Atom_such_c      = float(a[m1.start('Atom_such_c'):m1.end('Atom_such_c')])
            
            atom_record = lcccsObjects.Atom_Record()
            atom_record.set(ff)
            res = force_field.Add_Atom_Record(atom_record)
            
            print "load element", ff.Atom_element
            print ff.Atom_ff_type, ff.Atom_element, ff.Atom_such_n
#            atom_record.show_info()     


