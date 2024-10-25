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
 This test is to check that libint2, molint and GWP integrals are consistent
"""

import os
import sys
import math
import copy
import pytest
#import unittest

from liblibra_core import *

set1 = [  [1.0, 1.0,    0.0, 0.0, 0.0 ], 
          [1.0, 1.0,    1.0, 0.0, 0.0 ],
          [1.0, 1.0,    0.0, 1.0, 0.0 ],
          [1.0, 1.0,    0.0, 0.0, 1.0 ],
          [1.0, 1.0,    10.0, 0.0, 0.0 ],
          [1.0, 1.0,    0.0, 10.0, 0.0 ],
          [1.0, 1.0,    0.0, 0.0, 10.0 ],
          [0.1, 1.0,    1.0, 0.0, 0.0 ],
          [0.1, 1.0,    0.0, 1.0, 0.0 ],
          [0.1, 1.0,    0.0, 0.0, 1.0 ],
          [1.0, 0.1,    1.0, 0.0, 0.0 ],
          [1.0, 0.1,    0.0, 1.0, 0.0 ],
          [1.0, 0.1,    0.0, 0.0, 1.0 ]
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
    s = gaussian_overlap(0, 0, 0, a1, o1,  0, 0, 0, a2, o2)

    assert abs( x[0].get(0,0) - s ) < 1e-10



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
    #mu = transition_dipole_moment(0, 0, 0, a1, o1,  0, 0, 0, a2, o1, 1)
    sx = gaussian_overlap_ref(0, a1, 0.0, 0, a2, rx)
    sy = gaussian_overlap_ref(0, a1, 0.0, 0, a2, ry)
    sz = gaussian_overlap_ref(0, a1, 0.0, 0, a2, rz)
    
    nrm = 1.0/math.sqrt( gaussian_overlap_ref(0, a1, 0.0, 0, a1, 0.0)  * gaussian_overlap_ref(0, a2, 0.0, 0, a2, 0.0) )
    nrm = nrm**3    

    mux = gaussian_moment_ref(1, 0.0, 0.0,  0, a1, 0.0,  0, a2, rx) * sy * sz * nrm
    muy = gaussian_moment_ref(1, 0.0, 0.0,  0, a1, 0.0,  0, a2, ry) * sx * sz * nrm 
    muz = gaussian_moment_ref(1, 0.0, 0.0,  0, a1, 0.0,  0, a2, rz) * sx * sy * nrm

    assert abs( x[1].get(0,0) - mux ) < 1e-10 and abs( x[2].get(0,0) - muy ) < 1e-10 and abs( x[3].get(0,0) - muz ) < 1e-10


@pytest.mark.parametrize(('a1', 'a2', 'rx', 'ry', 'rz'), set1)
def test_pp(a1, a2, rx, ry, rz):

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
    #<g_a | [C0 + C2*(r-R)^2]*exp(-alp*(r-R)^2) | g_b>
    o = VECTOR(0.0, 0.0, 0.0)
    pp = pseudopot02(0.0, 1.0, 0.0, o,   0, 0, 0, a1, o1,  0, 0, 0, a2, o2)

    assert abs( x[4].get(0,0) + x[7].get(0,0) + x[9].get(0,0) - pp ) < 1e-10















