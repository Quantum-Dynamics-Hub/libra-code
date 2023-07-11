#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file vesta2qe.py 
# This module performs a conversion of a VESTA format files to the QE format. It makes sure 
# the fractional occupation atoms and other possibly overlapping atoms are removed (a unique one is left).


import re
import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def exclude(lab, R, min_dist, a,b,c,alp,bet,gam, out_file):
    # http://www.mse.mtu.edu/~drjohn/my3200/stereo/sg4.html

    c_alp = math.cos(math.radians(alp))
    s_alp = math.sin(math.radians(alp))

    c_bet = math.cos(math.radians(bet))
    s_bet = math.sin(math.radians(bet))

    c_gam = math.cos(math.radians(gam))
    s_gam = math.sin(math.radians(gam))


    M = MATRIX(3,3)
    M.set(0,0, a*s_bet);   M.set(0,1, b*s_alp*c_gam);  M.set(0,2, 0.0);
    M.set(1,0, 0.0);       M.set(1,1, b*s_alp*s_gam);  M.set(1,2, 0.0);
    M.set(2,0, a*c_bet);   M.set(2,1, b*c_alp);        M.set(2,2, c);

   

    M2 = M.T() * M

    sz = len(R)
    exclude_list = []

    for i in range(0,sz):

        if i not in exclude_list:
        
            for j in range(i+1, sz):

                for t1 in [-1.0, 0.0, 1.0]:
                    for t2 in [-1.0, 0.0, 1.0]:
                        for t3 in [-1.0, 0.0, 1.0]:

                            rij = R[i] - R[j]
                            rij.add(0, 0, t1)
                            rij.add(1, 0, t2)
                            rij.add(2, 0, t3)
            
                            d2 = (rij.T() * M2 * rij).get(0)
                            d = math.sqrt(d2)

                            if d < min_dist:
                                exclude_list.append(j)

#    print exclude_list

    types = []
    for i in range(0,sz):
        if i not in exclude_list:
            if lab[i] not in types:
                types.append(lab[i])

    


    #=========== Print the QE file =============
    f = open(out_file, "w")


    t1 = VECTOR(1.0, 0.0, 0.0)
    t2 = VECTOR(0.0, 1.0, 0.0)
    t3 = VECTOR(0.0, 0.0, 1.0)

    f.write("&SYSTEM\n")
    f.write("nat= %5i\n" % (len(R)-len(exclude_list)) )
    f.write("ntyp= %5i\n" % (len(types)))
    f.write("ibrav= 0\n" )
    f.write("/\n\n" )


    T1 = M * t1; T2 = M * t2; T3 = M * t3
    f.write("CELL_PARAMETERS Angstrom\n")
    f.write("%18.10f%16.10f%16.10f\n" % (T1.x, T1.y, T1.z))
    f.write("%18.10f%16.10f%16.10f\n" % (T2.x, T2.y, T2.z))
    f.write("%18.10f%16.10f%16.10f\n" % (T3.x, T3.y, T3.z))
    f.write("\nATOMIC_POSITIONS Angstrom\n")

    for i in range(0,sz):
        if i not in exclude_list:
            r = M * R[i]
            f.write("%2s%16.10f%16.10f%16.10f\n" % (lab[i], r.get(0), r.get(1), r.get(2)))

    f.close()




def parse_vesta(in_file, out_file, min_dist):

    """ -------- Load molecular system -------------
    Specify format manually (this gives flexibility) of format recognizable by program
    \param[in] out_file The name of the file containing the QChem output   

    """

    ######### Assume coordinates are given in Angstroms #############

    Angst_to_Bohr = 1.889725989  


    #------- Here are some basic patterns -------------
    INT    = '([1-9]([0-9]*))'
    NINT   = '([0-9]+)'
    SP     = '\s+'    
    DOUBLE = '([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
    WORD   = '([a-zA-Z]+)'
    ID     = '(([a-zA-Z]+)([a-zA-Z]+|\d+)*)'
    ID2     = '(([a-zA-Z]+|\d+)([a-zA-Z]+)*)'
    PHRASE = '"((\w|\W)+)"'
    compINT = re.compile(INT)


    #------- Here we define a format of file ----------
    # p - means 'Pattern'
    f   = open(in_file,'r')
    A = f.readlines()    
    f.close()
    sz = len(A)


    #============== Read in the cell parameters ========================
    pCell_keyword = '(?P<CELL_P>'+'CELLP'+')'+SP

    pcell_a = '(?P<cell_a>'+DOUBLE+')'+SP
    pcell_b = '(?P<cell_b>'+DOUBLE+')'+SP
    pcell_c = '(?P<cell_c>'+DOUBLE+')'+SP
    pcell_alp = '(?P<cell_alp>'+DOUBLE+')'+SP
    pcell_bet = '(?P<cell_bet>'+DOUBLE+')'+SP
    pcell_gam = '(?P<cell_gam>'+DOUBLE+')'+SP

    pCell_rec = pcell_a + pcell_b + pcell_c + pcell_alp + pcell_bet + pcell_gam

    a,b,c, alp, bet, gam = 0.0, 0.0, 0.0,  0.0, 0.0, 0.0

    for i in range(0,sz): 
        m1 = re.search(pCell_keyword,A[i])

        if m1!=None: 
            m2 = re.search(pCell_rec,A[i+1])

            if m2!=None:

                a = float(A[i+1][m2.start('cell_a'):m2.end('cell_a')])
                b = float(A[i+1][m2.start('cell_b'):m2.end('cell_b')])
                c = float(A[i+1][m2.start('cell_c'):m2.end('cell_c')])
                alp = float(A[i+1][m2.start('cell_alp'):m2.end('cell_alp')])
                bet = float(A[i+1][m2.start('cell_bet'):m2.end('cell_bet')])
                gam = float(A[i+1][m2.start('cell_gam'):m2.end('cell_gam')])


#    print "Cell parameters"
#    print a, b, c, alp, bet, gam


    #============== Read in the coordinates ========================
    pSTRUC_keyword = '(?P<STRUC>'+'STRUC'+')'+SP
    pAtom_id  = '(?P<Atom_id>'+WORD+')'+SP   
    pAtom_occ = '(?P<Atom_occ>'+DOUBLE+')'+SP
    pAtom_x = '(?P<Atom_x>'+DOUBLE+')'+SP
    pAtom_y = '(?P<Atom_y>'+DOUBLE+')'+SP
    pAtom_z = '(?P<Atom_z>'+DOUBLE+')'+SP

#    pAtom_rec = INT + pAtom_id + WORD + pAtom_occ + pAtom_x + pAtom_y + pAtom_z + ID2 + INT
    pAtom_rec = pAtom_id+ pAtom_occ + pAtom_x + pAtom_y + pAtom_z

    at_label, at_occ, at_coord = [], [], []

    for i in range(0,sz): 
        m1 = re.search(pSTRUC_keyword,A[i])

        if m1!=None: 
            nat, keep_going = 0, 1

            while keep_going:
                m2 = re.search(pAtom_rec,A[i+1+2*nat])

                if m2!=None:

                    at_label.append(  A[i+1+2*nat][m2.start('Atom_id'):m2.end('Atom_id')] )
                    at_occ.append(    float(A[i+1+2*nat][m2.start('Atom_occ'):m2.end('Atom_occ')]) )

                    x = float(A[i+1+2*nat][m2.start('Atom_x'):m2.end('Atom_x')])
                    y = float(A[i+1+2*nat][m2.start('Atom_y'):m2.end('Atom_y')])
                    z = float(A[i+1+2*nat][m2.start('Atom_z'):m2.end('Atom_z')])
                    r = MATRIX(3,1)
                    r.set(0, x); r.set(1, y); r.set(2, z)
                    at_coord.append( r )

                    nat = nat + 1
                else:
                    keep_going = 0


#    print "Atoms"
#    for i in xrange(len(at_label)):   
#        print i, at_label[i], at_occ[i], at_coord[i].show_matrix()


    exclude(at_label, at_coord, min_dist, a,b,c,alp,bet,gam, out_file)


# Example of usage
if __name__=='__main__':
    parse_vesta("1.vesta", "1.in", 0.7)
