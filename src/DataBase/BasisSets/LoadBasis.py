import lcccsObjects
import re
import os

def Parse(pattern,src):
# This function collect as many patterns 'pattern'
# in source string 'src' and the dictionary of gound groups is
# appended to 'out' array (list)
    out = []
    inp = src
    stat = 1
    while stat==1:
        m1 = re.search(pattern,inp)
        if m1!=None:
            out.append(m1.groupdict())           
            inp = inp[m1.end():]
            stat = 1
        else:
            stat = 0

    return out


def Generate_Orbital(Elt,type,coeffs,exponents):
# This function creates a single AO of given orbital moment, which is
# given by argument 'type' in terms of its expansion in linear combination
# of gaussian primitives. Expansion coefficients are given in argument
# 'coeffs' and the exponents of gaussian primitives are given in argument
# 'exponents'. This function return a STO AO object. Parameter 

    Orb_moments = {'S':0,'P':1,'D':2,'F':3}
    class aux:
        pass

    sto = lcccsObjects.AO()

    s = aux()
    s.element = Elt
    s.ao_shell = type
    sto.set(s)

    sz = len(exponents)
    i = 0
    while i<sz:
        g = aux()
        g.x_exp = Orb_moments[type]
        g.y_exp = 0
        g.z_exp = 0
        g.G_alpha = float(exponents[i])
        
        
        G = lcccsObjects.PrimitiveG()
        G.set(g)

        sto.add_primitive(float(coeffs[i]),G)

        i = i + 1

    sto.show_info()
    
    return sto





def Generate_Shell_Orbitals(sto):
# Generate all possivle orbitals with orbital moment = L
# from one such orbital
# That is given 2px orbital we can generate whole 2p shell = [2px,2py,2pz]
# Expansion of the radial part of all orbitals in the shell is the same
# and is given by expansion in argument 'sto', so we just copy that information

    class aux:
        pass

    Shell_orbitals = []

    L = sto.x_exp + sto.y_exp + sto.z_exp

    l = 0
    while l<=L:
        m = 0
        while m<=(L-l):
            n = 0
            while n<=(L-l-m):
                if l+m+n==L:
                    sto1 = lcccsObjects.AO()
                    sto1 = sto
                    sto1.x_exp = l
                    sto1.y_exp = m
                    sto1.z_exp = n

                    Shell_orbitals.append(sto1)

                n = n + 1
            m = m + 1
        l = l + 1

    return Shell_orbitals



def LoadBasis(basis_name):

    ff_file = os.getcwd() + "/DataBase/BasisSets/"+basis_name+".dat" # In current directory
    if os.path.isfile(ff_file):
        print "Loading basis from file = ",ff_file
    else:
        print "File ", ff_file, " does not exist"
        exit(1)
    

#------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'
    #DOUBLE = '([-+]?(\d*\.\d*)([eE][-+]?\d+)?)'
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    PHRASE = '"((\w|\W)+)"'
    NL     = "\\n"

#------- Here we define a format of file ----------
# p - means 'Pattern'
    pElt_symbol       = '(?P<Elt_symbol>'+WORD + ')' + SP
    pShell_type1      = '(?P<Shell_type1>S|P|D)'+SP 
    pShell_type2      = '(?P<Shell_type2>SP)'+ SP
    pG_exponent       = '(?P<G_exponent>'+DOUBLE+')'+SP
    pG_coeff1         = '(?P<G_coeff1>'+DOUBLE+')'+SP
    pG_coeff2         = '(?P<G_coeff2>'+DOUBLE+')'+SP


    Shell1_header=  pElt_symbol + pShell_type1
    Shell2_header=  pElt_symbol + pShell_type2
    Shell1_entry =  pG_exponent + pG_coeff1
    Shell2_entry =  pG_exponent + pG_coeff1 + pG_coeff2

    Shell1      = '(?P<Shell1>'+ Shell1_header + '('+Shell1_entry + SP +')+)'
    Shell2      = '(?P<Shell2>'+ Shell2_header + '('+Shell2_entry + SP +')+)'



    # Read data file
    f   = open(ff_file,'r')
    A = f.readlines()
    f.close() 
        

    # Concatenate all lines of the data file
    file_lines = ""
    for a in A:
        file_lines = file_lines + a
    a = file_lines

    class tmp:
        pass


    All_Shells1 = Parse(Shell1,a)
    All_Shells2 = Parse(Shell2,a)

    sz1 = len(All_Shells1)
    sz2 = len(All_Shells2)

#  Extracting S, P and D shells
    i = 0    
    while i<sz1:

        b = All_Shells1[i]['Shell1']
        print b

        Shell_h = Parse(Shell1_header,b)
        print Shell_h[0]
        Shell = Parse(Shell1_entry,b)
      
#----------------------------------------------------------------
        szx = len(Shell)
        A = []
        B = []

        j = 0
        while j<szx:
            A.append( float(Shell[j]['G_exponent']) )
            B.append( float(Shell[j]['G_coeff1']) )
            j = j + 1

        sto  = Generate_Orbital(Shell_h[0]['Elt_symbol'],Shell_h[0]['Shell_type1'],B,A)
        stos = Generate_Shell_Orbitals(sto)

        for stoi in stos:
            stoi.show_info()

#-----------------------------------------------------------------

        i = i + 1

# Extracting SP shells
    i = 0
    while i<sz2:
        b = All_Shells2[i]['Shell2']
        print b

        Shell_h = Parse(Shell2_header,b)
        print Shell_h[0]
        Shell = Parse(Shell2_entry,b)

        szx = len(Shell)
        A = []
        B = []
        C = []

        j= 0
        while j<szx:
            C.append( float(Shell[j]['G_coeff2']) )
            B.append( float(Shell[j]['G_coeff1']) )
            A.append( float(Shell[j]['G_exponent']) )            
            j = j + 1

        sto1  = Generate_Orbital(Shell_h[0]['Elt_symbol'],'S',B,A)
        stos = Generate_Shell_Orbitals(sto1)

        sto2  = Generate_Orbital(Shell_h[0]['Elt_symbol'],'P',C,A)
        stop = Generate_Shell_Orbitals(sto2)

        print len(stos)
        print len(stop)

#        for stoi in stos:
#            stoi.show_info()

        k = 0
        while k<len(stop):
            stop[k].show_info()
            k = k + 1


        i = i + 1







#LoadBasis("STO-3G")




