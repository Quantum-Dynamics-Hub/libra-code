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


set1 = [[1.0, 0.1],
        [10.0, 1.0],
        [0.0, 0.0],
        [-1.0, -0.1]
        ]


@pytest.mark.parametrize(('y', 'val'), set1)
def test_kcrpmd(y, val):

    ham = nHamiltonian(2, 2, 1)
    res = ham.kcrpmd_effective_potential(Py2Cpp_double([y]))

    assert abs(res - val) < 1e-10
