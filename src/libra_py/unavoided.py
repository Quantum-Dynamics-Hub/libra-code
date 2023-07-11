#*********************************************************************************  
#* Copyright (C) 2017 Kosuke Sato, Alexey V. Akimov 
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version. 
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>. 
#* 
#*********************************************************************************/
"""@package unavoided
 This module implements the functions to correct the state crossings (trivial, or 
 "unavoided") crossings.
 The order of the eigenstates is determined by the energies, not the orbital's character.         
 Say, if states i and j are localized on B and A, respectively, at "t+dt" while they were
 on A and B at "t", it is likely that this is just a manifestation of the trivial crossing, 
 not the actual transition. So one must correct the identity of the states.
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


    

def get_reordering(time_overlap):
    """ This function identifies which states have changed their identities via unavoided 
    (a.k.a. trivial) crossing when the system evolved in the interval from t to t+dt
    We start by looking at the time overlap matrix: <phi_i(t)|phi_i(t+dt)> 
    If no spurious state changes have happened during this time, the diagonal elements
    should be close to 1.0. If they are not - we locate to which state the transitions might 
    have happened.

    \param[in] time_overlap ( MATRIX ) the time overlap matrix, <phi_i(t)|phi_j(t+dt)>.
    Returns:
    perm - list of integers that describe the permutation. That is:
    perm[i] - is the index identifying the "older" state "i". Now, it may be labeled
    by some other index, j. 
    iperm[j] - if we were looking at a specific state labeled "j", iperm[j] this will be
    the index of that state in the present list of states (e.g. energies, etc.)

    """

    print("ERROR in libra_py.unavoided.get_reordering(...)")
    print("The implementation contained a bug. It is now re-implemented in C++\
          with the bug fixed. Please use the C++ version. Exiting...")
    sys.exit(0)
    

    # extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    S = CMATRIX(time_overlap)   # just a temporary working object
    sz = time_overlap.num_of_rows
    perm = range(sz)  # original permutation
    

    for col in range(0,sz):

        indx = -1
        val = 0.0+0.0j
        while indx!=col:

            # Find the max element in the given column "col"
            [indx, val] = S.max_col_elt(col)
            
            # Apply permutation (col, indx) to the present "perm" list
            tmp = perm[col]
            perm[col] = perm[indx]
            perm[indx] = tmp

            # Do the corresponding swap of the columns in the S matrix
            S.swap_cols(col,indx)

    return perm



def _test_setup():
    """Prepare the test matrices"""

    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j); 
    a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 
    #print "a = ";  a.show_matrix()

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j); 
    b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);
    #print "b = ";  b.show_matrix()


    # non-identical (4x4) matrix
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 1 0 0 0 |
    # | 0 0 0 1 |
    c = CMATRIX(4,4)
    c.set(0,1,1.0+0j);  c.set(1,2,1.0+0j); 
    c.set(2,0,1.0+0j);  c.set(3,3,1.0+0j);
    #print "b = ";  b.show_matrix()


    # non-identical (8x8) matrix
    # | 1 0 0 0 0 0 0 0 |
    # | 0 1 0 0 0 0 0 0 |
    # | 0 0 0 1 0 0 0 0 |
    # | 0 0 1 0 0 0 0 0 |
    # | 0 0 0 0 1 0 0 0 |
    # | 0 0 0 0 0 0 0 1 |
    # | 0 0 0 0 0 1 0 0 |
    # | 0 0 0 0 0 0 1 0 |

    d = CMATRIX(8,8)
    d.set(0,0,1.0+0j);  d.set(1,1,1.0+0j);
    d.set(2,3,1.0+0j);  d.set(3,2,1.0+0j);
    d.set(4,4,1.0+0j);  d.set(5,7,1.0+0j);
    d.set(6,5,1.0+0j);  d.set(7,6,1.0+0j);

    #print "c = ";   c.show_matrix()

    return a, b, c


class TestUnavoided(unittest.TestCase):
    def test_reordering(self):
        """Tests the reordering algorithm"""
        a,b,c = _test_setup()

        perm_a = get_reordering(a)
        print("Input matrix "); a.show_matrix()
        print("Permutation = ", perm_a)
        self.assertEqual(perm_a, [0,1,2,3])

        perm_b = get_reordering(b)
        print("Input matrix "); b.show_matrix()
        print("Permutation = ", perm_b)
        self.assertEqual(perm_b, [0,2,3,1])

        perm_c = get_reordering(c)
        print("Input matrix "); c.show_matrix()
        print("Permutation = ", perm_c)
        self.assertEqual(perm_c, [1,2,0,3])


if __name__=='__main__':
    unittest.main()

