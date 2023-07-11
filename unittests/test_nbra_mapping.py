#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math
import time

from liblibra_core import *
import libra_py.workflows.nbra.mapping as mapping


S = CMATRIX(6,6) 
#     1                           2                      3                    4                      5                     6
#  H-1, a                         H,a                   L, a                 H-1, b                  H, b                 L, b
S.set(0,0, 1.0+0.0j);   S.set(0, 1, 2.0+0.0j);  S.set(0,2, 3.0+0.0j);  S.set(0,3, 4.0+0.0j);  S.set(0,4, 5.0+0.0j);  S.set(0,5, 6.0+0.0j);   # 1 
S.set(1,0, 2.0+0.0j);   S.set(1, 1,12.0+0.0j);  S.set(1,2,13.0+0.0j);  S.set(1,3,14.0+0.0j);  S.set(1,4,15.0+0.0j);  S.set(1,5,16.0+0.0j);   # 2
S.set(2,0, 3.0+0.0j);   S.set(2, 1,13.0+0.0j);  S.set(2,2,23.0+0.0j);  S.set(2,3,24.0+0.0j);  S.set(2,4,25.0+0.0j);  S.set(2,5,26.0+0.0j);   # 3
S.set(3,0, 4.0+0.0j);   S.set(3, 1,14.0+0.0j);  S.set(3,2,24.0+0.0j);  S.set(3,3,34.0+0.0j);  S.set(3,4,35.0+0.0j);  S.set(3,5,36.0+0.0j);   # 4
S.set(4,0, 5.0+0.0j);   S.set(4, 1,15.0+0.0j);  S.set(4,2,25.0+0.0j);  S.set(4,3,35.0+0.0j);  S.set(4,4,45.0+0.0j);  S.set(4,5,46.0+0.0j);   # 5
S.set(5,0, 6.0+0.0j);   S.set(5, 1,16.0+0.0j);  S.set(5,2,26.0+0.0j);  S.set(5,3,36.0+0.0j);  S.set(5,4,46.0+0.0j);  S.set(5,5,56.0+0.0j);   # 6

#==================== Test 1 ===========================
print("Test 1")
# Simple case of single excitation of the same spin
GS = [1, -1, 2, -2]
S1 = [1, -1, 3, -2]
print(mapping.ovlp_arb_mo(GS, S1, S) )  # ok - pyxaid style

GSa = [1, -4, 2, -5]
S1a = [1, -4, 3, -5]
print(mapping.ovlp_arb_mo(GSa, S1a, S) )  


#==================== Test 2 ===========================
print("Test 2")
# Simple case of single excitation of the same spin, alpha for now
GS = [1, -1, 2, -2]
S1 = [1, -1, 2, -3]
print(mapping.ovlp_arb_mo(GS, S1, S) )  # ok - pyxaid style

# Simple case of single excitation of the same spin, analogously, but in the beta channel
GSa = [1, -4, 2, -5]
S1a = [1, -4, 2, -6]
print(mapping.ovlp_arb_mo(GSa, S1a, S) )  


#==================== Test 3 ===========================
print("Test 3")
# Coupling of the spin-flipped configurations
S1 = [1, -1, 2, -3]
S2 = [1, -1, 2,  3]
print(mapping.ovlp_arb_mo(S1, S2, S) )  # ok - pyxaid style

S1 = [1, -4, 2, -6]
S2 = [1, -4, 2,  3]
print(mapping.ovlp_arb_mo(S1, S2, S) )  # ok - pyxaid style


#==================== Test 4 ===========================
print("Test 4")
# Double excitation
GS = [1, -1, 2, -2]
D1 = [3, -1, 2, -3]
print(mapping.ovlp_arb_mo(GS, D1, S) )  # ok - pyxaid style






