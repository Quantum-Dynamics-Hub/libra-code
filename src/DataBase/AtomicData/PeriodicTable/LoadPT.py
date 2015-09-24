import lcccsObjects
import re

def Load_PT(ObjSpc,pt_file):
       
#    pt_file = "./DataBase/AtomicData/PeriodicTable/elements.dat" # The path relative to the root (where lcccsObjects.py sits)

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
    pElement_name                = '(?P<Element_name>'+WORD+')'+SP
    pElement_atomic_number       = '(?P<Element_atomic_number>'+INT+')'+SP
    pElement_nucleus_charge      = '(?P<Element_nucleus_charge>'+INT+')'+SP
    pElement_number_of_electrons = '(?P<Element_number_of_electrons>'+INT+')'+SP
    pElement_atomic_mass         = '(?P<Element_atomic_mass>'+DOUBLE+')'+SP
    pElement_period              = '(?P<Element_period>'+INT+')'+SP
    pElement_group               = '(?P<Element_group>'+INT+')'+SP
    pElement_block               = '(?P<Element_block>'+WORD+')'+SP
    pElement_color               = '(?P<Element_color>'+INT+')'+SP
    pElement_bond_length         = '(?P<Element_bond_length>'+DOUBLE+')'+SP
    pElement_atomic_radius       = '(?P<Element_atomic_radius>'+DOUBLE+')'+SP


    Element_Record = pElement_name + \
                 pElement_atomic_number + \
                 pElement_nucleus_charge + \
                 pElement_number_of_electrons + \
                 pElement_atomic_mass + \
                 pElement_period + \
                 pElement_group + \
                 pElement_block + \
                 pElement_color + \
                 pElement_bond_length + \
                 pElement_atomic_radius



#-------- Load periodic system --------------------
# Defines atom type by element name

    f   = open(pt_file,'r')
    A = f.readlines()
    f.close()

    class element:
        pass

    for a in A:
        m1 = re.search(Element_Record,a)
        if m1!=None:
            elem = element()
            elem.Elt_name        = a[m1.start('Element_name'):m1.end('Element_name')]
            elem.Elt_id          = int(float(a[m1.start('Element_atomic_number'):m1.end('Element_atomic_number')]))
            elem.Elt_number      = elem.Elt_id
            elem.Elt_nucleus_charge = elem.Elt_id
            elem.Elt_mass        = float(a[m1.start('Element_atomic_mass'):m1.end('Element_atomic_mass')])
            elem.Elt_period      = int(float(a[m1.start('Element_period'):m1.end('Element_period')]))
            elem.Elt_group       = int(float(a[m1.start('Element_group'):m1.end('Element_group')]))
            elem.Elt_block       = a[m1.start('Element_block'):m1.end('Element_block')]
#            elem.Atom_color       = int(float(a[m1.start('Element_color'):m1.end('Element_color')]))
            # put bond_length here
            elem.Elt_atomic_radius = float(a[m1.start('Element_atomic_radius'):m1.end('Element_atomic_radius')])


           # print "load  element", element.Atom_element

            elt_record = lcccsObjects.Element()
            elt_record.set(elem)
            res = ObjSpc.Add_Element_To_Periodic_Table(elt_record)
 
            print "load element", elem.Elt_name
            elt_record.show_info()
          #  print res


