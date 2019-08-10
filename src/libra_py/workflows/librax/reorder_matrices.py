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
## \file reorder_matrices.py 
# This module implements a function that reorders elements of density matrices 
# or energy ones according to a permutation "perm". 
# Also it includes unittest on the function with some examples.

import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
import unavoided_tmp # This will be deleted when it is merged into unavoided.py of libra_py

def reorder(p,A,E):
    # This function reorders columns of density matrix "A" and diagonal elements of energy "E"
    # following permutation "p".
    # Note: reordering columns of "A" means reordering newer eigenvectors "sd_basis2". 
    #
    # param[in]     p a list of permutation indices
    # param[in,out] A MATRIX object including density matrix <phi_i(t)|phi_j(t+dt)>
    # param[in,out] E MATRIX object including energies in diagonal elements

    perm = p[:] # create a working list

    for col in xrange(len(perm)):
        indx=perm[col]
        if indx!=col:
            A.swap_cols(indx,col)
            E.swap_cols(indx,col); E.swap_rows(indx,col)

            indx2 = perm[indx]
            if indx2!=col:
                perm[col] = col
                perm[indx2] = indx
            else:
                perm[col] = col
                perm[indx] = indx
            #print "col is %i and indx is %i" % (col,indx)
            #print "temporal perm_c is"; print perm_c
            #print "temporal matrix is"; c.show_matrix()

def _prepare_density_matrices():
    
    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j);  a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j);  b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);

    # non-identical (6x6) matrix
    # | 0 0 1 0 0 0|
    # | 0 0 0 1 0 0|
    # | 1 0 0 0 0 0|
    # | 0 0 0 0 0 1|
    # | 0 0 0 0 1 0|
    # | 0 1 0 0 0 0|
    c = CMATRIX(6,6)
    c.set(0,2,1.0+0j);  c.set(1,3,1.0+0j);  c.set(2,0,1.0+0j);  c.set(3,5,1.0+0j);
    c.set(4,4,1.0+0j);  c.set(5,1,1.0+0j);

    # non-identical (7x7) matrix (This contains 3 groups of doubly mixed states)
    # | 1    0    0    0    0    0    0 |
    # | 0 0.71 0.72    0    0    0    0 |
    # | 0 0.73 0.74    0    0    0    0 |
    # | 0    0    0 0.72 0.73    0    0 |
    # | 0    0    0 0.74 0.75    0    0 |
    # | 0    0    0    0    0 0.76 0.77 |
    # | 0    0    0    0    0 0.78 0.79 |
    d = CMATRIX(7,7)
    d.set(0,0,1.0+0j);
    d.set(1,1,0.71+0j); d.set(1,2,0.72+0j);
    d.set(2,1,0.73+0j); d.set(2,2,0.74+0j); 
    d.set(3,3,0.72+0j); d.set(3,4,0.73+0j);
    d.set(4,3,0.74+0j); d.set(4,4,0.75+0j);
    d.set(5,5,0.76+0j); d.set(5,6,0.77+0j);
    d.set(6,5,0.78+0j); d.set(6,6,0.79+0j);

    return a,b,c,d

def _prepare_energy_matrices():

    # | 0  0  0  0 |
    # | 0 11  0  0 |
    # | 0  0 22  0 |
    # | 0  0  0 33 |
    Ea = CMATRIX(4,4)
    Ea.set(0,0,0.0+0j); Ea.set(1,1,11.0+0j); Ea.set(2,2,22.0+0j); Ea.set(3,3,33.0+0j);
    
    Eb = CMATRIX(Ea)

    # | 0  0  0  0  0  0 |
    # | 0 11  0  0  0  0 |
    # | 0  0 22  0  0  0 |
    # | 0  0  0 33  0  0 |
    # | 0  0  0  0 44  0 |
    # | 0  0  0  0  0 55 |

    Ec = CMATRIX(6,6) # eigenenergy matrix (diagonal)
    Ec.set(0,0,0.0+0j); Ec.set(1,1,11.0+0j); Ec.set(2,2,22.0+0j); Ec.set(3,3,33.0+0j);
    Ec.set(4,4,44.0+0j); Ec.set(5,5,55.0+0j);

    # | 0  0  0  0  0  0  0 |
    # | 0 11  0  0  0  0  0 |
    # | 0  0 22  0  0  0  0 |
    # | 0  0  0 33  0  0  0 |
    # | 0  0  0  0 44  0  0 |
    # | 0  0  0  0  0 55  0 |
    # | 0  0  0  0  0  0 66 |

    Ed = CMATRIX(7,7) # eigenenergy matrix (diagonal)
    Ed.set(0,0,0.0+0j); Ed.set(1,1,11.0+0j); Ed.set(2,2,22.0+0j); Ed.set(3,3,33.0+0j);
    Ed.set(4,4,44.0+0j); Ed.set(5,5,55.0+0j); Ed.set(6,6,66.0+0j);

    return Ea,Eb,Ec,Ed

class Test_unavoided(unittest.TestCase):
    def test_reorder(self):

        rnd = Random()

        a,b,c,d = _prepare_density_matrices()
        Ea,Eb,Ec,Ed = _prepare_energy_matrices()
        #sys.exit(0)

        p4 = range(a.num_of_rows) # [0,1,2,3]
        p6 = range(c.num_of_rows) # [0,1,2,3,4,5]
        p7 = range(d.num_of_rows) # [0,1,2,3,4,5,6]

        print "p4 is ",p4
        print "p6 is ",p6
        print "p7 is ",p7

        '''extract indices for reordering '''
        #perm_a = unavoided.get_reordering(a)
        perm_a = unavoided_tmp.get_reordering(a)
        print "Input density matrix a"; a.show_matrix()
        print "Permutation a = ", perm_a
        self.assertEqual(perm_a, [0,1,2,3])

        #perm_b = unavoided.get_reordering(b)
        perm_b = unavoided_tmp.get_reordering(b)
        print "Input density matrix b"; b.show_matrix()
        print "Input energy matrix Eb" ; Eb.show_matrix()
        print "Permutation b = ", perm_b
        self.assertEqual(perm_b, [0,2,3,1])

        #perm_c = unavoided.get_reordering(c)
        perm_c = unavoided_tmp.get_reordering(c)
        print "Input density matrix c"; c.show_matrix()
        print "Input energy matrix Ec"; Ec.show_matrix()
        print "Permutation c = ", perm_c
        self.assertEqual(perm_c, [2,3,0,5,4,1])

        perm_d = unavoided_tmp.get_reordering(d)
        print "Input density matrix d"; d.show_matrix()
        print "Input energy matrix Ed"; Ed.show_matrix()
        print "Permutation d = ", perm_d
        perm_d_eq = [0,1,2,3,4,5,6,0,1,2,3,4,6,5,0,1,2,4,3,5,6,0,1,2,4,3,6,5,0,2,1,3,4,5,6,0,2,1,3,4,6,5,0,2,1,4,3,5,6,0,2,1,4,3,6,5]
        self.assertEqual(perm_d, perm_d_eq)

        #sys.exit(0)
        ''' permutation according to "perm" list '''

        if p4 != perm_a:
            reorder(perm_a,a,Ea); print "Matrix a is reordered"
        print "Output density matrix a"; a.show_matrix()
        print "Output energy matrix Ea"; Ea.show_matrix()
        print "Permutation a = ", perm_a

        if p4 != perm_b:
            reorder(perm_b,b,Eb); print "Matrix b is reordered"
        print "Output density matrix b"; b.show_matrix()
        print "Output energy matrix Eb"; Eb.show_matrix()
        print "Permutation b = ", perm_b

        if p6 != perm_c:
            reorder(perm_c,c,Ec); print "Matrix c is reordered"
        print "Output density matrix c"; c.show_matrix()
        print "Output energy matrix Ec"; Ec.show_matrix()
        print "Permutation c = ", perm_c

        if p7 != perm_d:
            ksi = rnd.uniform(0.0,1.0)
            Np = len(perm_d)/d.num_of_rows # the number of permutations
            ip = int(Np*ksi)
            istart = ip*d.num_of_rows
            iend = (ip+1)*d.num_of_rows
            ptmp = perm_d[istart:iend]
            reorder(ptmp,d,Ed); print "Matrix d is reordered"
            print "ip is %i" % ip
        print "Output density matrix d"; d.show_matrix()
        print "Output energy matrix Ed"; Ed.show_matrix()
        print "Permutation d = ", perm_d


if __name__=='__main__':
    unittest.main()
    
