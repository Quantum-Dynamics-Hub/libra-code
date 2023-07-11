#*********************************************************************************
#* Copyright (C) 2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
#
# This example demonstrates how to get normal modes info from QE (phonon) output 
# 
###################################################################

import os
import math
import sys
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import QE_methods as qem


Elts, R, U = qem.get_QE_normal_modes("silicon.dyn1",1)
Elts, R, U = qem.get_QE_normal_modes("Cs4SnBr6_T200.dyn1")

