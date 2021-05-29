#*********************************************************************************
#* Copyright (C) 2017-2019 Brendan Smith, Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: autoconnect
   :platform: Unix, Windows
   :synopsis: This module implements functions to determine the connectivity 
       matrix in molecular systems.

.. moduleauthor:: Brendan Smith, Ekadashi Pradhan, Alexey V. Akimov

"""


import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
import util.libutil as comn

def autoconnect(R, MaxCoord, params):
    """

    Args:
        R ( list of VECTOR objects ): The atomic coordinates of the system [ units: arbitrary ]
        MaxCoord ( list of ints ): Maximal coorination numbers of each atom 
        params ( dictionary ): Parameters controlling the execution of this function

            * **params["Rcut"]** ( double ): The maximal radius of connectivity [ units: same as R, default: 0.0 ]
            * **params["tv1"]** ( VECTOR ): unit cell translation vector a [ units: same as R ]
            * **params["tv2"]** ( VECTOR ): unit cell translation vector b [ units: same as R ]
            * **params["tv3"]** ( VECTOR ): unit cell translation vector c [ units: same as R ]
            * **params["pbc_opt"]** ( string ): what type of periodicity to assume. 
                Options: "a", "b", "c", "ab", "ac", "bc", "abc", "none" [ default: None ]
            * **params["opt"]** ( int ): Option to obey maximal coordination number: 
                - 0: the connected atoms in a pair have to obey both the coordination numbers
                   Some of the atoms may stay undercoordinated. This is the consistent scheme in
                   the sense that if A is connected to B, then B is necessarily connected to A [default]
                - 1: treat the MaxCoord as the minimal coordination number each atom should attain.
                   Some atoms may be overcoordinated. 
            * **params["verbosity"]** ( int ): Flag to control the amount of the printout [ default: 0 ]

    Returns:
        tuple: ( res, line, unsorted_pairs_out ), where:

            * res ( list of lists [ res[i][0], res[0][1] ] ), where 
                * res[i][0] ( int ): index of the atom
                * res[i][1] ( list of ints ): indices of the atoms that are connected to res[i][0]

            * line ( string ): the connectivity information in the .ent format

            * unsorted_pairs_out ( same type as ```res```): same as ```res```, but not sorted
                                  
   
    """

    critical_params = [ ] 
    default_params = { "Rcut":0.0, "pbc_opt":"none", "opt":0,
                       "tv1":VECTOR(1.0e+10, 0.0, 0.0),
                       "tv2":VECTOR(0.0, 1.0e+10, 0.0),
                       "tv3":VECTOR(0.0, 0.0, 1.0e+10),
                       "verbosity":0
                     }
    comn.check_input(params, default_params, critical_params)

    Rcut = params["Rcut"]
    pbc_opt = params["pbc_opt"]
    verbosity = params["verbosity"]
    opt = params["opt"]
    tv1 = params["tv1"]
    tv2 = params["tv2"]
    tv3 = params["tv3"]


    # Preprocessing - sanity check
    N = len(R)

    if len(MaxCoord) != N:
        print("Error: The length of the MaxCoord list should be the same as the lenght of the R list!")
        print("Exiting now...")
        sys.exit(0)

 
    unsorted_pairs = []
    mapping = []
    periodicity = []    


    if pbc_opt not in ["a", "b", "c", "ab", "ac", "bc", "abc","none"]:
        print("Error: pbc_opt ", pbc_opt, " is not recognized")
        sys.exit(0)
    
    transl_a = [0]
    transl_b = [0]
    transl_c = [0]

    if pbc_opt in ["a", "ab", "ac", "abc"]:
        transl_a = [-1, 0, 1]
    if pbc_opt in ["b", "ab", "bc", "abc"]:
        transl_b = [-1, 0, 1]
    if pbc_opt in ["c", "ac", "bc", "abc"]:
        transl_c = [-1, 0, 1]


 
    # Distances between all the pairs
    count = 0
    for n1 in transl_a:
        for n2 in transl_b:
            for n3 in transl_c:

                for i in range(0,N):    

                    start = 0
                    if(n1==0 and n2==0 and n3==0):
                        start = i+1

                    for j in range(start,N):        

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
        print("Pairs of atoms separated no more than by Rcut")
        print(unsorted_pairs)
        print("Formatted printout:")
        for it in unsorted_pairs:
            print(F"Atoms {mapping[it[0]][0]} and {mapping[it[0]][1]} are separated by\
                  {it[1]:8.5f}. Periodic translation applied = {periodicity[it[0]][0]:8.5f},\
                  {periodicity[it[0]][1]:8.5f}, {periodicity[it[0]][2]:8.5f} )" )
        print("Pairs of atoms separated no more than by Rcut, sorted by the distance")
        print(sorted_pairs)
        print("Formatted printout:")
        for it in sorted_pairs:
            print(F"Atoms {mapping[it[0]][0]} and {mapping[it[0]][1]} are separated by\
                  {it[1]:8.5f}. Periodic translation applied = {periodicity[it[0]][0]:8.5f},\
                  {periodicity[it[0]][1]:8.5f}, {periodicity[it[0]][2]:8.5f} )" )



    # Initialize the results variables
    act_coord = [0]*N   # actual coordination numbers
    num_sat = 0         # the number of valence-saturated atoms
    res = []
    for i in range(0,N):
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
    for i in range(0,N):
        res[i][1] = sorted(res[i][1])
   

    line = ""
    for it in res:
        line = line + "CONECT %5i " % (it[0]+1)
        for it2 in it[1]:
            line = line + " %5i " % (it2+1)
        line = line + "\n"


    if verbosity==1:
        print("res = ", res)
        print("Formatted output:")
        print(line)

 
    return res, line, unsorted_pairs_out



def find_undercoordinated_atoms(res, MaxCoord):
    """

    Args:
        res ( list of lists [ res[i][0], res[0][1] ] ), where 
            * res[i][0] ( int ): index of the atom
            * res[i][1] ( list of ints ): indices of the atoms that are connected to res[i][0]

            This would typically be the first output of the ```autoconnect``` function

        MaxCoord ( list of ints ): Maximal coorination numbers of each atom 

    Returns: 
        (list):  out, such that:

            out[i][0] - is the index of the i-th undercoordinated atom
            out[i][1] - is the number of dangling bonds on the i-th atom

    """

    nat = len(MaxCoord)  # the total number of atoms
    out = []

    for i in range(0,nat):
        diff = MaxCoord[i] - len(res[i][1])
        if  diff > 0:
            out.append([i, diff])

    return out


def example_1():

    rnd = Random()
    N = 10
    R = []
    MaxCoord = []

    for i in range(0,N):
        if rnd.uniform(0.0, 1.0)<0.5:
            MaxCoord.append(3)
        else:
            MaxCoord.append(3)
        R.append(VECTOR(rnd.normal(), rnd.normal(), rnd.normal()))
    Rcut = 3.0

    print("Runing   autoconnect(R, MaxCoord, Rcut, 0, 1)  with : ")
    print("R = ", R)
    print("MaxCoord = ", MaxCoord)
    print("Rcut = ", Rcut)
    
    autoconnect(R, MaxCoord, {"Rcut":Rcut, "opt":0, "verbosity":1})



class TestAutoconnect(unittest.TestCase):
    """ Docs
    """

    def test_1(self):
        """ Diatomic"""

        R = []
        MaxCoord = []
        R.append(VECTOR(0.0, 0.0, 0.0));   MaxCoord.append(1)  # Atom 1
        R.append(VECTOR(1.0, 0.0, 0.0));   MaxCoord.append(1)  # Atom 2

                 
        res, lines = autoconnect(R, MaxCoord, {"Rcut":1.5})          
        self.assertEqual(res, [[0, [1]], [1, [0]]] )

        res, lines = autoconnect(R, MaxCoord, {"Rcut":0.5})          
        self.assertEqual(res, [[0,[]], [1, []]] )


    def test_2(self):
        """ Triatomic """

        R = []
        MaxCoord = []
        R.append(VECTOR(0.0, 0.0, 0.0));   MaxCoord.append(2)  # Atom 1
        R.append(VECTOR(1.0, 0.0, 0.0));   MaxCoord.append(2)  # Atom 2
        R.append(VECTOR(1.2, 0.0, 0.0));   MaxCoord.append(2)  # Atom 2

                 
        res, lines = autoconnect(R, MaxCoord, {"Rcut":1.5})          
        self.assertEqual(res, [[0, [1,2]], [1, [0,2]], [2, [0,1]] ] )

        res, lines = autoconnect(R, MaxCoord, {"Rcut":0.5})          
        print(res)
        self.assertEqual(res, [ [0, []], [1, [2]], [2, [1]] ] )





if __name__=='__main__':    
    unittest.main()
    #example_1()
