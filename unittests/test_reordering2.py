#*********************************************************************************  
#* Copyright (C) 2017-2018 Kosuke Sato, Alexey V. Akimov 
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version. 
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>. 
#* 
#*********************************************************************************/
"""
 This is an extended test suite for validating the reordering algorithm implemented in C++
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

    # non-identical matrix (corresponding to doubly mixed states)
    # | 1    0    0 0 |
    # | 0 0.76 0.79 0 |
    # | 0 0.79 0.80 0 |
    # | 0    0    0 1 |

    e = CMATRIX(4,4)
    e.set(0,0,1.0+0j);
    e.set(1,1,0.76+0j); e.set(1,2,0.79+0j);
    e.set(2,1,0.79+0j); e.set(2,2,0.80+0j);
    e.set(3,3,1.0+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 |
    # | 0 0.71 0.75 0.77 0 |
    # | 0 0.72 0.75 0.78 0 |
    # | 0 0.73 0.76 0.80 0 |
    # | 0    0    0    0 1 |

    f = CMATRIX(5,5)
    f.set(0,0,1.0+0j);
    f.set(1,1,0.71+0j); f.set(1,2,0.75+0j); f.set(1,3,0.77+0j);
    f.set(2,1,0.72+0j); f.set(2,2,0.75+0j); f.set(2,3,0.78+0j);
    f.set(3,1,0.73+0j); f.set(3,2,0.76+0j); f.set(3,3,0.80+0j);
    f.set(4,4,1.0+0j);

    # non-identical matrix (corresponding to 3 groups of doubly mixed states)
    # | 1    0    0    0    0    0    0|
    # | 0 0.73 0.75    0    0    0    0|
    # | 0 0.76 0.76    0    0    0    0|
    # | 0    0    0 0.71 0.72    0    0|
    # | 0    0    0 0.74 0.78    0    0|
    # | 0    0    0    0    0 0.77 0.78|
    # | 0    0    0    0    0 0.79 0.81|
    
    g = CMATRIX(7,7)
    g.set(0,0,1.0+0j);
    g.set(1,1,0.73+0j); g.set(1,2,0.75+0j);
    g.set(2,1,0.76+0j); g.set(2,2,0.76+0j);
    g.set(3,3,0.71+0j); g.set(3,4,0.72+0j);
    g.set(4,3,0.74+0j); g.set(4,4,0.78+0j);
    g.set(5,5,0.77+0j); g.set(5,6,0.78+0j);
    g.set(6,5,0.79+0j); g.set(6,6,0.81+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 |
    # | 0 0.77 0.78 0.78 0 |
    # | 0 0.76 0.75 0.89 0 |
    # | 0 0.75 0.73 0.77 0 |
    # | 0    0    0    0 1 | 

    h = CMATRIX(5,5)
    h.set(0,0,1.0+0j);
    h.set(1,1,0.77+0j); h.set(1,2,0.78+0j); h.set(1,3,0.78+0j);
    h.set(2,1,0.76+0j); h.set(2,2,0.75+0j); h.set(2,3,0.89+0j);
    h.set(3,1,0.75+0j); h.set(3,2,0.73+0j); h.set(3,3,0.77+0j);
    h.set(4,4,1.00+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 | 
    # | 0 0.65 0.78 0.77 0 |
    # | 0 0.89 0.75 0.76 0 |
    # | 0 0.77 0.54 0.65 0 |
    # | 0    0    0    0 1 |

    h1 = CMATRIX(5,5)
    h1.set(0,0,1.0+0j);
    h1.set(1,3,0.77+0j); h1.set(1,2,0.78+0j); h1.set(1,1,0.65+0j);
    h1.set(2,3,0.76+0j); h1.set(2,2,0.75+0j); h1.set(2,1,0.89+0j);
    h1.set(3,3,0.65+0j); h1.set(3,2,0.54+0j); h1.set(3,1,0.77+0j);
    h1.set(4,4,1.00+0j);

    return a, b, c, d, e, f, g, h, h1


class TestUnavoided(unittest.TestCase):
    def test_reordering(self):
        """Tests the reordering algorithm"""
        a,b,c,d,e,f,g,h, h1 = _test_setup()

        perm_a = Cpp2Py(get_reordering(a))
        print "Input matrix "; a.show_matrix()
        print "Permutation = ", perm_a
        self.assertEqual(perm_a, [0,1,2,3])

        perm_b = Cpp2Py(get_reordering(b))
        print "Input matrix "; b.show_matrix()
        print "Permutation = ", perm_b
#        self.assertEqual(perm_b, [0,2,3,1])
        self.assertEqual(perm_b, [0,3,1,2])

        perm_c = Cpp2Py(get_reordering(c))
        print "Input matrix "; c.show_matrix()
        print "Permutation = ", perm_c
        self.assertEqual(perm_c, [1,2,0,3])

        perm_e = get_reordering(e)
        print "Input matrix "; e.show_matrix()
        print "Permutation = ", perm_e
        self.assertEqual(perm_e, [0,2,1,3])

        perm_f = get_reordering(f)
        print "Input matrix "; f.show_matrix()
        print "Permutation = ", perm_f
        perm_f_eq = [0, 2, 1, 3, 4]
        self.assertEqual(perm_f, perm_f_eq)

        perm_g = get_reordering(g)
        print "Input matrix "; g.show_matrix()
        print "Permutation = ", perm_g
        perm_g_eq = [0, 2, 1, 3, 4, 5, 6]
        self.assertEqual(perm_g, perm_g_eq)

        perm_h = get_reordering(h)
        print "Input matrix "; h.show_matrix()
        print "Permutation = ", perm_h
        perm_h_eq = [0, 2, 3, 1, 4]
        self.assertEqual(perm_h, perm_h_eq)

        perm_h1 = get_reordering(h1)
        print "Input matrix "; h1.show_matrix()
        print "Permutation = ", perm_h1
        perm_h1_eq = [0, 2, 1, 3, 4]
        self.assertEqual(perm_h1, perm_h1_eq)

if __name__=='__main__':
    unittest.main()
