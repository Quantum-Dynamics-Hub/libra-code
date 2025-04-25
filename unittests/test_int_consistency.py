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
 This test is to check that libint2, molint and GWP integrals are consistent
"""

import os
import sys
import math
import copy
import pytest
# import unittest

from liblibra_core import *

set1 = [[1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, 1.0],
        [1.0, 1.0, 10.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 10.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, 10.0],
        [0.1, 1.0, 1.0, 0.0, 0.0],
        [0.1, 1.0, 0.0, 1.0, 0.0],
        [0.1, 1.0, 0.0, 0.0, 1.0],
        [1.0, 0.1, 1.0, 0.0, 0.0],
        [1.0, 0.1, 0.0, 1.0, 0.0],
        [1.0, 0.1, 0.0, 0.0, 1.0]
        ]


@pytest.mark.parametrize(('a1', 'a2', 'rx', 'ry', 'rz'), set1)
def test_ovlps(a1, a2, rx, ry, rz):

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(rx, ry, rz)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    # Format of x:
    # ['S', 'x', 'y', 'z', 'x2', 'xy', 'xz', 'y2', 'yz', 'z2', 'x3', 'x2y', 'x2z', 'xy2', 'xyz', 'xz2', 'y3', 'y2z', 'yz2', 'z3']
    x = compute_emultipole3(s1, s2, 1)

    # Internal molints:
    s = gaussian_overlap(0, 0, 0, a1, o1, 0, 0, 0, a2, o2)

    assert abs(x[0].get(0, 0) - s) < 1e-10


@pytest.mark.parametrize(('a1', 'a2', 'rx', 'ry', 'rz'), set1)
def test_dipoles(a1, a2, rx, ry, rz):

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(rx, ry, rz)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    # Format of x:
    # ['S', 'x', 'y', 'z', 'x2', 'xy', 'xz', 'y2', 'yz', 'z2', 'x3', 'x2y', 'x2z', 'xy2', 'xyz', 'xz2', 'y3', 'y2z', 'yz2', 'z3']
    x = compute_emultipole3(s1, s2, 1)

    # Internal molints:
    mu = transition_dipole_moment(0, 0, 0, a1, o1, 0, 0, 0, a2, o2, 1)

    assert abs(x[1].get(0, 0) - mu.x) < 1e-10 and abs(x[2].get(0, 0) -
                                                      mu.y) < 1e-10 and abs(x[3].get(0, 0) - mu.z) < 1e-10


@pytest.mark.parametrize(('a1', 'a2', 'rx', 'ry', 'rz'), set1)
def test_dipoles2(a1, a2, rx, ry, rz):
    """
    More explicit calculation of the dipole moments
    """

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(rx, ry, rz)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    # Format of x:
    # ['S', 'x', 'y', 'z', 'x2', 'xy', 'xz', 'y2', 'yz', 'z2', 'x3', 'x2y', 'x2z', 'xy2', 'xyz', 'xz2', 'y3', 'y2z', 'yz2', 'z3']
    x = compute_emultipole3(s1, s2, 1)

    # Internal molints:
    sx = gaussian_overlap(0, a1, 0.0, 0, a2, rx, 1)
    sy = gaussian_overlap(0, a1, 0.0, 0, a2, ry, 1)
    sz = gaussian_overlap(0, a1, 0.0, 0, a2, rz, 1)

    mux = gaussian_moment(0, a1, 0.0, 1, 0.0, 0.0, 0, a2, rx, 1) * sy * sz
    muy = gaussian_moment(0, a1, 0.0, 1, 0.0, 0.0, 0, a2, ry, 1) * sx * sz
    muz = gaussian_moment(0, a1, 0.0, 1, 0.0, 0.0, 0, a2, rz, 1) * sx * sy

    assert abs(x[1].get(0, 0) - mux) < 1e-10 and abs(x[2].get(0, 0) -
                                                     muy) < 1e-10 and abs(x[3].get(0, 0) - muz) < 1e-10


@pytest.mark.parametrize(('a1', 'a2', 'rx', 'ry', 'rz'), set1)
def test_quadratic(a1, a2, rx, ry, rz):

    e1 = Py2Cpp_double([a1])
    c1 = Py2Cpp_double([1.0])
    o1 = VECTOR(0.0, 0.0, 0.0)
    s1 = initialize_shell(0, 1, e1, c1, o1)

    e2 = Py2Cpp_double([a2])
    c2 = Py2Cpp_double([1.0])
    o2 = VECTOR(rx, ry, rz)
    s2 = initialize_shell(0, 1, e2, c2, o2)

    # Format of x:
    # ['S', 'x', 'y', 'z', 'x2', 'xy', 'xz', 'y2', 'yz', 'z2', 'x3', 'x2y', 'x2z', 'xy2', 'xyz', 'xz2', 'y3', 'y2z', 'yz2', 'z3']
    x = compute_emultipole3(s1, s2, 1)

    # Internal molints:
    # Internal molints:
    sx = gaussian_overlap(0, a1, 0.0, 0, a2, rx, 1)
    sy = gaussian_overlap(0, a1, 0.0, 0, a2, ry, 1)
    sz = gaussian_overlap(0, a1, 0.0, 0, a2, rz, 1)

    x2 = gaussian_moment(0, a1, 0.0, 2, 0.0, 0.0, 0, a2, rx, 1) * sy * sz
    y2 = gaussian_moment(0, a1, 0.0, 2, 0.0, 0.0, 0, a2, ry, 1) * sx * sz
    z2 = gaussian_moment(0, a1, 0.0, 2, 0.0, 0.0, 0, a2, rz, 1) * sx * sy

    x3 = gaussian_moment(0, a1, 0.0, 3, 0.0, 0.0, 0, a2, rx, 1) * sy * sz
    y3 = gaussian_moment(0, a1, 0.0, 3, 0.0, 0.0, 0, a2, ry, 1) * sx * sz
    z3 = gaussian_moment(0, a1, 0.0, 3, 0.0, 0.0, 0, a2, rz, 1) * sx * sy

    assert abs(x[4].get(0, 0) -
               x2) < 1e-10 and abs(x[7].get(0, 0) -
                                   y2) < 1e-10 and abs(x[9].get(0, 0) -
                                                       z2) < 1e-10 and abs(x[10].get(0, 0) -
                                                                           x3) < 1e-10 and abs(x[16].get(0, 0) -
                                                                                               y3) < 1e-10 and abs(x[19].get(0, 0) -
                                                                                                                   z3) < 1e-10
