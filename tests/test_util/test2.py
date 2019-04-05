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

import os
import math
import sys

#if sys.platform=="cygwin":
#    from cyglibra_core import *
#elif sys.platform=="linux" or sys.platform=="linux2":
#    from liblibra_core import utils.libutils

#from libra_py import *
import util.libutil as comn


params = {"a":1.0, "b":"A" }
default_params = {"c":"23d" }

print params
comn.check_input(params, default_params, [] )
print params

comn.check_input(params, default_params, ["d"] )
