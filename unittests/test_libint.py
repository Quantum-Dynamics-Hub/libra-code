# *********************************************************************************
# * Copyright (C) 2024 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
 This is an extended test suite for validating the reordering algorithm implemented in C++
"""

import os
import sys
import math
import copy
import pytest
# import unittest

from liblibra_core import *


def test_initialize_shell():
    # Signature:
    # initialize_shell(int l_val, bool is_spherical, const
    # std::vector<double>& exponents, const std::vector<double>& coeff,
    # VECTOR& coords){

    exps = Py2Cpp_double([1.0])
    coeff = Py2Cpp_double([1.0])
    coord = VECTOR(0.0, 0.0, 0.0)

    s1 = initialize_shell(0, 1, exps, coeff, coord)

    assert len(s1) == 1


set1 = [[1.0, 1.0, 0.0, 1.0]]


@pytest.mark.parametrize(('a1', 'a2', 'r', 'val'), set1)
def test_overlaps(a1, a2, r, val):

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(r, 0.0, 0.0)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    x = compute_overlaps(s1, s2, 1)

    assert abs(x.get(0, 0) - val) < 1e-10


set2 = [[1.0, 1.0, 0.0, [1.0, 0.0, 0.0, 0.0]]]


@pytest.mark.parametrize(('a1', 'a2', 'r', 'vals'), set2)
def test_multipoles(a1, a2, r, vals):

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(r, 0.0, 0.0)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    x = compute_emultipole3(s1, s2, 1)

    res = True
    for i in range(4):
        res = res * (abs(x[i].get(0, 0) - vals[i]) < 1e-10)
    assert res


# ['S', 'x', 'y', 'z', 'x2', 'xy', 'xz', 'y2', 'yz', 'z2', 'x3', 'x2y', 'x2z', 'xy2', 'xyz', 'xz2', 'y3', 'y2z', 'yz2', 'z3']
