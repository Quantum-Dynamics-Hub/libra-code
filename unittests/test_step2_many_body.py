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
import libra_py.workflows.nbra.step2_many_body as mb


set1 = [(68, [67, 68, 69, 70, 71, 72], [[[68, 69], "alp"], [[68, 69], "bet"]], 1, 0, [[1, 2, -7, -8], [1, 3, -7, -8], [1, 2, -7, -9]]),
        (68, [67, 68, 69, 70, 71, 72], [[[68, 69], "alp"], [[68, 69], "bet"]], 2, 0, [[1, -1, 2, -2], [1, -1, 3, -2], [1, -1, 2, -3]])
        ]


@pytest.mark.parametrize(('ks_orbital_homo_index', 'ks_orbital_indicies', 'sd_basis_states',
                         'sd_format', 'ks_beta_homo_index', 'sd_basis'), set1)
def test_reindex_cp2k_sd_states(
        ks_orbital_homo_index,
        ks_orbital_indicies,
        sd_basis_states,
        sd_format,
        ks_beta_homo_index,
        sd_basis):

    res = mb.reindex_cp2k_sd_states(
        ks_orbital_homo_index,
        ks_orbital_indicies,
        sd_basis_states,
        sd_format,
        ks_beta_homo_index)

    for i in range(len(sd_basis)):
        assert res[i] == sd_basis[i]
