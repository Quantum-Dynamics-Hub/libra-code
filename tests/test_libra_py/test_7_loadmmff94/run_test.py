#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
# Tutorial: This file demonstrates the functionality of the 
# "LoadTRIPOS" module of the "libra_py"
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *



mmff94 = ForceField()

LoadMMFF94.Load_MMFF94(mmff94,"mmff94/mmff94_types1.dat", "mmff94/mmff94_types2.dat",
                 "mmff94/mmff94_bonds1.dat", "mmff94/mmff94_bonds2.dat","mmff94/mmff94_bonds3.dat",
                 "mmff94/mmff94_angles1.dat","mmff94/mmff94_angles2.dat","mmff94/mmff94_angles3.dat",
                 "mmff94/mmff94_torsions.dat", "mmff94/mmff94_oop.dat" )


mmff94.show_info()
