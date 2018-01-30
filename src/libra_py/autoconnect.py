#*********************************************************************************
#* Copyright (C) 2017 Brendan Smith, Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  The code that automatically determines a connectivity matrix in an arbitrary molecular system.
#
#

import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#from libra_py import *



def autoconnect(R, MaxCoord, Rcut, opt=0, verbosity=0):
    """
    \param[in] R (list of VECTORs) The atomic coordinates of the system
    \param[in] MaxCoord (list of ints) The list of maxima coorination numbers
    \param[in] Rcut (list of floats) Radius of inclusion for each atom 
    \param[in] opt (int) Option to obey maximal coordination number: 
               0 - (default) the connected atoms in a pair have to obey both coordination numbers
                   Some of the atoms may stay undercoordinated. This is the consistent scheme in
                   the sense that if A is connected to B, then B is necessarily connected to A
               1 - treat the MaxCoord as the minimal coordination number each atom should attain.
                   Some atoms may be overcoordinated. 

    \param[in] verbosity (int) Flag to control the amount of the printout
    Returns:
   
    """


    # Preprocessing - sanity check
    N = len(R)

    if len(MaxCoord) != N:
        print "Error: The length of the MaxCoord list should be the same as the lenght of the R list!"
        print "Exiting now..."
        sys.exit(0)



 
    unsorted_pairs = []
    mapping = []
    

    # Distances between all the pairs
    count = 0
    for i in range(0,N):    
        for j in range(i+1,N):        
            r = (R[i]-R[j]).length()

            if r<Rcut:
                unsorted_pairs.append([count, r])
                mapping.append([i,j])
                count += 1

    # Sort all the pairs according to the interparticle distance
    sorted_pairs = merge_sort(unsorted_pairs) 


    if verbosity==1:
        print "Pairs of atoms separated no more than by Rcut"
        print unsorted_pairs
        print "Formatted printout:"
        for it in unsorted_pairs:
            print "Atoms %5i and %5i are separated by %8.5f " % (mapping[it[0]][0], mapping[it[0]][1], it[1])

        print "Pairs of atoms separated no more than by Rcut, sorted by the distance"
        print sorted_pairs
        print "Formatted printout:"
        for it in sorted_pairs:
            print "Atoms %5i and %5i are separated by %8.5f " % (mapping[it[0]][0], mapping[it[0]][1], it[1])



    # Initialize the results variables
    act_coord = [0]*N   # actual coordination numbers
    num_sat = 0         # the number of valence-saturated atoms
    res = []
    for i in xrange(N):
        res.append([i,[]])

    # Now lets pick the pairs with the minimal distance, depending on the option...
    for it in sorted_pairs:
        if num_sat < N:        
            i,j = mapping[it[0]][0], mapping[it[0]][1]

            if opt==0:  #... no overcoordinated atoms
                if act_coord[i]<MaxCoord[i] and act_coord[j]<MaxCoord[j]:  
                    res[i][1].append(j)
                    res[j][1].append(i)
                    act_coord[i] += 1
                    act_coord[j] += 1

                if act_coord[i]==MaxCoord[i]:
                    num_sat += 1

                if act_coord[j]==MaxCoord[j]:
                    num_sat += 1



            elif opt==1: # ... no undercoordinated atoms
                if act_coord[i]<MaxCoord[i]:  
                    res[i][1].append(j)
                    act_coord[i] += 1

                if act_coord[j]<MaxCoord[j]:  
                    res[j][1].append(i)
                    act_coord[j] += 1

    # Order the indices of the atoms to which each atom is connected (for testing and convenience)    
    for i in xrange(N):
        res[i][1] = sorted(res[i][1])
   

    line = ""
    for it in res:
        line = line + "CONECT %5i " % (it[0]+1)
        for it2 in it[1]:
            line = line + " %5i " % (it2+1)
        line = line + "\n"


    if verbosity==1:
        print "res = " , res
        print "Formatted output:"
        print line

 
    return res, line





def autoconnect_pbc(R, MaxCoord, Rcut, tv1, tv2, tv3, pbc_opt, opt=0, verbosity=0):
    """
    \param[in] R (list of VECTORs) The atomic coordinates of the system
    \param[in] MaxCoord (list of ints) The list of maxima coorination numbers
    \param[in] Rcut (list of floats) Radius of inclusion for each atom 
    \param[in] tv1 (VECTOR) - unit cell translation vector a
    \param[in] tv2 (VECTOR) - unit cell translation vector b
    \param[in] tv3 (VECTOR) - unit cell translation vector c
    \param[in] pbc_opt (string) - "a", "b", "c", "ab", "ac", "bc", "abc"
    \param[in] opt (int) Option to obey maximal coordination number: 
               0 - (default) the connected atoms in a pair have to obey both coordination numbers
                   Some of the atoms may stay undercoordinated. This is the consistent scheme in
                   the sense that if A is connected to B, then B is necessarily connected to A
               1 - treat the MaxCoord as the minimal coordination number each atom should attain.
                   Some atoms may be overcoordinated. 

    \param[in] verbosity (int) Flag to control the amount of the printout
    Returns:
   
    """


    # Preprocessing - sanity check
    N = len(R)

    if len(MaxCoord) != N:
        print "Error: The length of the MaxCoord list should be the same as the lenght of the R list!"
        print "Exiting now..."
        sys.exit(0)

 
    unsorted_pairs = []
    mapping = []
    periodicity = []    


    if pbc_opt not in ["a", "b", "c", "ab", "ac", "bc", "abc"]:
        print "Error: pbc_opt ", pbc_opt, " is not recognized"
        sys.exit(0)
    
    transl_a = [0.0]
    transl_b = [0.0]
    transl_c = [0.0]

    if pbc_opt in ["a", "ab", "ac", "abc"]:
        transl_a = [-1.0, 0.0, 1.0]
    if pbc_opt in ["b", "ab", "bc", "abc"]:
        transl_b = [-1.0, 0.0, 1.0]
    if pbc_opt in ["c", "ac", "bc", "abc"]:
        transl_c = [-1.0, 0.0, 1.0]


 
    # Distances between all the pairs
    count = 0
    for i in range(0,N):    
        for j in range(i+1,N):        
            for n1 in transl_a:
                for n2 in transl_b:
                    for n3 in transl_c:

                        T = n1 * tv1 + n2 * tv2 + n3 * tv3
                        r = (R[i]-R[j]-T).length()

                        if r<Rcut:
                            unsorted_pairs.append([count, r])
                            mapping.append([i,j])
                            periodicity.append([n1,n2,n3])
                            count += 1

    # Sort all the pairs according to the interparticle distance
    sorted_pairs = merge_sort(unsorted_pairs) 


    unsorted_pairs_out = []
    for it in unsorted_pairs:
        unsorted_pairs_out.append( [ mapping[it[0]][0], mapping[it[0]][1],
                                   VECTOR(0.0, 0.0, 0.0),
                                   VECTOR( periodicity[it[0]][0] * tv1 + periodicity[it[0]][1] * tv2 + periodicity[it[0]][2] * tv3 ) ] 
                                  )
    
    if verbosity==1:
        print "Pairs of atoms separated no more than by Rcut"
        print unsorted_pairs
        print "Formatted printout:"
        for it in unsorted_pairs:
            print "Atoms %5i and %5i are separated by %8.5f. Periodic translation applied = [%8.5f, %8.5f, %8.5f] " \
            % (mapping[it[0]][0], mapping[it[0]][1], it[1], periodicity[it[0]][0], periodicity[it[0]][1], periodicity[it[0]][2] )

        print "Pairs of atoms separated no more than by Rcut, sorted by the distance"
        print sorted_pairs
        print "Formatted printout:"
        for it in sorted_pairs:
            print "Atoms %5i and %5i are separated by %8.5f. Periodic translation applied = [%8.5f, %8.5f, %8.5f] " \
            % (mapping[it[0]][0], mapping[it[0]][1], it[1], periodicity[it[0]][0], periodicity[it[0]][1], periodicity[it[0]][2] )



    # Initialize the results variables
    act_coord = [0]*N   # actual coordination numbers
    num_sat = 0         # the number of valence-saturated atoms
    res = []
    for i in xrange(N):
        res.append([i,[]])

    # Now lets pick the pairs with the minimal distance, depending on the option...
    for it in sorted_pairs:
        if num_sat < N:        
            i,j = mapping[it[0]][0], mapping[it[0]][1]

            if opt==0:  #... no overcoordinated atoms
                if act_coord[i]<MaxCoord[i] and act_coord[j]<MaxCoord[j]:  
                    res[i][1].append(j)
                    res[j][1].append(i)
                    act_coord[i] += 1
                    act_coord[j] += 1

                if act_coord[i]==MaxCoord[i]:
                    num_sat += 1

                if act_coord[j]==MaxCoord[j]:
                    num_sat += 1



            elif opt==1: # ... no undercoordinated atoms
                if act_coord[i]<MaxCoord[i]:  
                    res[i][1].append(j)
                    act_coord[i] += 1

                if act_coord[j]<MaxCoord[j]:  
                    res[j][1].append(i)
                    act_coord[j] += 1

    # Order the indices of the atoms to which each atom is connected (for testing and convenience)    
    for i in xrange(N):
        res[i][1] = sorted(res[i][1])
   

    line = ""
    for it in res:
        line = line + "CONECT %5i " % (it[0]+1)
        for it2 in it[1]:
            line = line + " %5i " % (it2+1)
        line = line + "\n"


    if verbosity==1:
        print "res = " , res
        print "Formatted output:"
        print line

 
    return res, line, unsorted_pairs_out



def example_1():

    rnd = Random()
    N = 10
    R = []
    MaxCoord = []

    for i in xrange(N):
        if rnd.uniform(0.0, 1.0)<0.5:
            MaxCoord.append(3)
        else:
            MaxCoord.append(3)
        R.append(VECTOR(rnd.normal(), rnd.normal(), rnd.normal()))
    Rcut = 3.0

    print "Runing   autoconnect(R, MaxCoord, Rcut, 0, 1)  with : "
    print "R = ", R
    print "MaxCoord = ", MaxCoord
    print "Rcut = ", Rcut
    
    autoconnect(R, MaxCoord, Rcut, 0, 1)



class TestAutoconnect(unittest.TestCase):
    """ Docs
    """

    def test_1(self):
        """ Diatomic"""

        R = []
        MaxCoord = []
        R.append(VECTOR(0.0, 0.0, 0.0));   MaxCoord.append(1)  # Atom 1
        R.append(VECTOR(1.0, 0.0, 0.0));   MaxCoord.append(1)  # Atom 2

                 
        res, lines = autoconnect(R, MaxCoord, 1.5)          
        self.assertEqual(res, [[0, [1]], [1, [0]]] )

        res, lines = autoconnect(R, MaxCoord, 0.5)          
        self.assertEqual(res, [[0,[]], [1, []]] )


    def test_2(self):
        """ Triatomic """

        R = []
        MaxCoord = []
        R.append(VECTOR(0.0, 0.0, 0.0));   MaxCoord.append(2)  # Atom 1
        R.append(VECTOR(1.0, 0.0, 0.0));   MaxCoord.append(2)  # Atom 2
        R.append(VECTOR(1.2, 0.0, 0.0));   MaxCoord.append(2)  # Atom 2

                 
        res, lines = autoconnect(R, MaxCoord, 1.5)          
        self.assertEqual(res, [[0, [1,2]], [1, [0,2]], [2, [0,1]] ] )

        res, lines = autoconnect(R, MaxCoord, 0.5)          
        print res
        self.assertEqual(res, [ [0, []], [1, [2]], [2, [1]] ] )





if __name__=='__main__':    
    unittest.main()
    #example_1()
