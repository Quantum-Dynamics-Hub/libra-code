# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
 This is module for testing mapping procedure for SD overlap calculations, etc.
"""

import os
import sys
import math
import copy
#import pytest
# import unittest
import pytest
import numpy as np

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

from libra_py.workflows.nbra.mapping import sd2indx

def generate():
    """
    Run this function to see what is the current results 
    or if you want to add more reference results
    """
    for inp in [ [1,-1, 2, -2], 
                 [-1,1, 2, -2],
                 [1,-1, 2, -3],
                 [1,-1,-2, 3]
               ]:
        for do_sort in [False, True]:
            for user_notation in [0, 1]:
                out, p = sd2indx(inp, nbasis=6, do_sort=do_sort, user_notation=user_notation)
                print(F"{inp} do_sort={do_sort}, user_notation={user_notation} -> {out}, {p}" )


# Test cases
test_data = [
    # inp, do_sort, user_notation, expected_out, expected_p
    ([1, -1, 2, -2], False, 0, [0, 0, 1, 1], 1),
    ([1, -1, 2, -2], False, 1, [0, 3, 1, 4], 1),
    ([1, -1, 2, -2], True, 0,  [0, 0, 1, 1], 1),
    ([1, -1, 2, -2], True, 1,  [0, 1, 3, 4], -1),

    ([-1, 1, 2, -2], False, 0, [0, 0, 1, 1], 1),
    ([-1, 1, 2, -2], False, 1, [3, 0, 1, 4], 1),
    ([-1, 1, 2, -2], True, 0,  [0, 0, 1, 1], 1),
    ([-1, 1, 2, -2], True, 1,  [0, 1, 3, 4], 1),

    ([1, -1, 2, -3], False, 0, [0, 0, 1, 2], 1),
    ([1, -1, 2, -3], False, 1, [0, 3, 1, 5], 1),
    ([1, -1, 2, -3], True, 0,  [0, 0, 1, 2], 1),
    ([1, -1, 2, -3], True, 1,  [0, 1, 3, 5], -1),

    ([1, -1, -2, 3], False, 0, [0, 0, 1, 2], 1),
    ([1, -1, -2, 3], False, 1, [0, 3, 4, 2], 1),
    ([1, -1, -2, 3], True,  0, [0, 0, 1, 2], 1),
    ([1, -1, -2, 3], True,  1, [0, 2, 3, 4], 1),
]


@pytest.mark.parametrize(
    "inp,do_sort,user_notation,expected_out,expected_p",
    test_data
)
def test_sd2indx(inp, do_sort, user_notation, expected_out, expected_p):
    """Test sd2indx mapping for several reference spin–determinant mappings."""
    out, p = sd2indx(inp, nbasis=6, do_sort=do_sort, user_notation=user_notation)

    # convert numpy array to list if needed
    if isinstance(out, np.ndarray):
        out = out.tolist()

    assert out == expected_out, f"Output mismatch for inp={inp}, do_sort={do_sort}, user_notation={user_notation}"
    assert p == expected_p, f"Phase mismatch for inp={inp}, do_sort={do_sort}, user_notation={user_notation}"

