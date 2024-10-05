#*********************************************************************************  
#* Copyright (C) 2017-2024 Kosuke Sato, Alexey V. Akimov 
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
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
import pytest
#import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def init_set():
    """Prepare the test matrices"""

    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j); 
    a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 
    a_perm = [0,1,2,3]
    #print "a = ";  a.show_matrix()

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j); 
    b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);
    b_perm = [0,2,3,1]
    #print "b = ";  b.show_matrix()

    # non-identical (4x4) matrix
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 1 0 0 0 |
    # | 0 0 0 1 |
    c = CMATRIX(4,4)
    c.set(0,1,1.0+0j);  c.set(1,2,1.0+0j); 
    c.set(2,0,1.0+0j);  c.set(3,3,1.0+0j);
    c_perm = [1,2,0,3]
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
    e_perm = [0,2,1,3]

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
    f_perm = [0,2,1,3,4]

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
    g_perm = [0,2,1,3,4,5,6]

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
    h_perm = [0,2,3,1,4]

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
    h1_perm = [0,2,1,3,4]

    return [ [a, a_perm], [b, b_perm], [c, c_perm], [e, e_perm], 
             [f, f_perm], [g, g_perm], [h, h_perm], [h1, h1_perm] ]

set1 = init_set()



@pytest.mark.parametrize(('mtx', 'perm'), set1)
def test_Munkres_Kunh(mtx, perm):
    en = CMATRIX(mtx); en *= 0.0;
    x = Cpp2Py( Munkres_Kuhn(mtx, en, 0.0, 0) )
    assert perm == x 

@pytest.mark.parametrize(('mtx', 'perm'), set1)
def test_get_reordering(mtx, perm):
    x = Cpp2Py( get_reordering(mtx) )
    assert perm == x

@pytest.mark.parametrize(('mtx', 'perm'), set1)
def test_hungarian(mtx, perm):
    en = CMATRIX(mtx); en *= 0.0;
    x = Cpp2Py( hungarian_algorithm(mtx, en, 0.0, 0) )
    assert perm == x



def get_U(N):
    U = CMATRIX(N,N)
    for i in range(N):
        U.set(i,i, 1.0+0.0j)
    return U

def permute(U, perm):
    pU = CMATRIX(U)
    pU.permute_cols(Py2Cpp_int(perm))
    return pU


set2 = [ [0,1], [1,0], 
         [0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,1,0], [2,0,1],
         [0,1,2,3], [0,2,1,3], [1,0,2,3], [1,2,0,3], [2,1,0,3], [2,0,1,3], 
         [1,2,3,0], [1,3,2,0], [2,1,3,0], [2,3,1,0], [3,2,1,0], [3,2,1,0],
         [3,1,2,0], [2,3,0,1], [2,0,3,1], [3,2,0,1], [3,0,2,1], [0,3,2,1],
         [0,2,3,1], [3,0,1,2], [3,1,0,2], [0,3,1,2], [0,1,3,2], [1,0,3,2],
         [1,3,0,2] 
       ]

@pytest.mark.parametrize(('perm'), set2)
def test_Munkres_Kunh_2(perm):
    N = len(perm)
    U  = get_U(N)
    Ut = permute(U, perm)
    X = U.H() * Ut
    en = CMATRIX(X); en *= 0.0;
    res = Cpp2Py( Munkres_Kuhn(X, en, 0.0, 0) )
    assert res == perm

@pytest.mark.parametrize(('perm'), set2)
def test_get_reordering_2(perm):
    N = len(perm)
    U  = get_U(N)
    Ut = permute(U, perm)
    X = Ut.H() * U
    res = Cpp2Py(get_reordering(X))
    assert res == perm


@pytest.mark.parametrize(('perm'), set2)
def test_hungarian(perm):
    N = len(perm)
    U  = get_U(N)
    Ut = permute(U, perm)
    X = Ut.H() * U
    en = CMATRIX(X); en *= 0.0;
    res = Cpp2Py(hungarian_algorithm(X, en, 0.0))
    assert res == perm


