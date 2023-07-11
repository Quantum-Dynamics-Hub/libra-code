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
# "init_system" module of the "libra_py": Lets try to construct multiple copies
###################################################################

import sys
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


rnd = Random()


label = ["H"]
R = [VECTOR()]
g = [VECTOR()]
T = 300.0
sigma = 0.1
df = 1

syst = []

f = open("system.xyz","w")
f.close()

# all excitations for each nuclear configuration
ninit = 5
for i in xrange(ninit):
    df = 0 # debug flag
    # Here we use libra_py module!
    # Utilize the gradients on the ground (0) excited state    
    syst.append( init_system.init_system(label, R, g, rnd, T, sigma, df, os.getcwd()+"/elements.txt") )    
    print syst[i]
    syst[i].print_xyz("system.xyz",i)


