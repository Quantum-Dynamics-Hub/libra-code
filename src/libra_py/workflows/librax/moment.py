#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file moment.py
# This program implements the module that calculates and returns
# the dipole moment matrixes at given space coordinates r like <MO(t)| r |MO(t+dt)>.

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def transition_dipole_moments(ao,C,act):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao : atomic orbital basis
    # \param[in] C  : MO-LCAO coefficients
    # \param[in] act : active_space
    #
    # Used in: gamess_to_libra.py/gamess_to_libra

    v = VECTOR(0.0,0.0,0.0) # 
    gx = PrimitiveG(1,0,0, 0.0, v) # = x
    gy = PrimitiveG(0,1,0, 0.0, v) # = y
    gz = PrimitiveG(0,0,1, 0.0, v) # = z
    g = [gx,gy,gz]

    Norb = len(ao)
    #mu = [MATRIX(Norb, Norb)]*3
    d = MATRIX(Norb,Norb)

    mu = [CMATRIX(Norb, Norb)]*3

    for k in xrange(3): # all components

        # moment matrices in the AO basis
        for i in xrange(Norb): # all orbitals
            for j in xrange(Norb):
                d.set(i,j,gaussian_moment(ao[i], g[k], ao[j]) )

        dd = CMATRIX(d)
        mu[k] = C.T() * dd * C
        #mu[k] = C.T() * d * C
    
    return mu[0], mu[1], mu[2]
