#*********************************************************************************
#* Copyright (C) 2024 Alexey V. Akimov
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

s1 = CMATRIX(2,2); s1.set(0,0, 1.0+0.0j); s1.set(1,1, 1.0+0.0j);
a1 = CMATRIX(2,2); a1.set(0,0, 1.0+0.0j); a1.set(1,1, 1.0+0.0j);
b1 = CMATRIX(2,2); b1.set(0,0, 1.0+0.0j); b1.set(1,1, 1.0+0.0j);

set1 = [ [ s1, a1, b1] ]

@pytest.mark.parametrize(('S', 'S_half', 'S_i_half'), set1)
def test_sqrt_matrix1(S, S_half, S_i_half):

    s_half = CMATRIX(S)
    s_i_half = CMATRIX(S)
    _id = CMATRIX(S); _id.identity()
    sqrt_matrix(S, s_half, s_i_half)
    #s_half.show_matrix()
    #s_i_half.show_matrix()

    assert abs((s_half - S_half).max_elt()) < 1e-10 and abs((s_i_half - S_i_half).max_elt()) <1e-10


set2 = [0.0, 45.0, 90.0]

@pytest.mark.parametrize(('phi'), set2)
def test_sqrt_matrix2(phi):
    """
    Some of these tests fail because the inverse matrices are defined up to a phase factor
    """

    phi = phi / 2.0*math.pi # in rad
    cs = math.cos(phi)*(1.0+0.0j)
    si = math.sin(phi)*(1.0+0.0j)
    A = CMATRIX(2,2);  A.set(0,0,cs); A.set(0,1,si); A.set(1,0,-si); A.set(1,1,cs)
    B = CMATRIX(2,2);  B.set(0,0,cs); B.set(0,1,-si); B.set(1,0,si); B.set(1,1,cs)
    S = A * A;
    s_i_half = CMATRIX(S)
    s_half = CMATRIX(S)
    sqrt_matrix(S, s_half, s_i_half)
    s_half.show_matrix()
    s_i_half.show_matrix()
    A.show_matrix()
    B.show_matrix()

    assert abs((A - s_half).max_elt()) < 1e-10 and abs((B - s_i_half).max_elt()) <1e-10


set2 = [0.0, 45.0, 90.0]

@pytest.mark.parametrize(('phi'), set2)
def test_orthogonalized_T(phi):
    """
    Despite the S^{1/2} and S^{-1/2} can be defined up to a complex phase,
    the Lowdin orthogonalization of the already orthogonal matrix should 
    give the starting matrix - so let's test that this is the case
    """
 
    phi = phi / 2.0*math.pi # in rad
    cs = math.cos(phi)*(1.0+0.0j)
    si = math.sin(phi)*(1.0+0.0j)
    T = CMATRIX(2,2);  T.set(0,0,cs); T.set(0,1,si); T.set(1,0,-si); T.set(1,1,cs)
    T2 = T.H() * T
    #B = CMATRIX(2,2);  B.set(0,0,cs); B.set(0,1,-si); B.set(1,0,si); B.set(1,1,cs)
    #S = A * A;
    T2_i_half = CMATRIX(T)
    T2_half = CMATRIX(T)
    sqrt_matrix(T2, T2_half, T2_i_half)
    T_ort = T * T2_i_half
    #s_half.show_matrix()
    #s_i_half.show_matrix()
    #A.show_matrix()
    #B.show_matrix()

    assert abs((T - T_ort).max_elt()) < 1e-10


