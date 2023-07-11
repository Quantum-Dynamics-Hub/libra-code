#***********************************************************
# * Copyright (C) 2018 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/


import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../workflows/pyxaid2")
sys.path.insert(1,cwd+"/../../workflows/pyxaid2")

import trajectory
import utils

#===================== Reading info ==================================
# Method 0 
# no SOC, non-polarized, one kpt
# To get HOMO, divide the number of electrons (72) by 2 => 36. Everything else is relative to it
print "================ no SOC, non-polarized, one kpt ==================="

params = {"nac_method":0, "prefix": "nosoc_nospinpol_1x1x1/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
E = e[0].show_matrix()

print "Active space orbitals"
print "uncomment"
coeff[0].show_matrix()

print "Verify the orthogonality of orbitals"
C = coeff[0]
S = C.H() * C
S.show_matrix()

print "Orthogonalize orbitals"
Ctilda = utils.orthogonalize_orbitals(C)
S = Ctilda.H() * Ctilda
S.show_matrix()




# no SOC, non-polarized, several kpts
print "================ no SOC, non-polarized, several kpts ==================="

params = {"nac_method":0, "prefix": "nosoc_nospinpol_2x2x2/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
sz = len(e)
for k in xrange(sz):
    print "==== kpt %i ======" % (k)

    print "Energies of the active space orbitals, Ha"
    E = e[k].show_matrix()

    print "Active space orbitals"
    print "uncomment"
    #coeff[k].show_matrix()

    print "Verify the orthogonality of orbitals"
    C = coeff[k]
    S = C.H() * C
    S.show_matrix()

    print "Orthogonalize orbitals"
    Ctilda = utils.orthogonalize_orbitals(C)
    S = Ctilda.H() * Ctilda
    S.show_matrix()




# Method 1
# no SOC, spin-polarized, one kpt
# To get HOMO, divide the number of electrons (72) by 2 => 36. Everything else is relative to it
print "================ no SOC, spin-polarized, one kpt ==================="

params = {"nac_method":1, "prefix": "nosoc_spinpol_1x1x1/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
E = e[0].show_matrix()
E = e[1].show_matrix()

print "Active space orbitals"
print "uncomment"
#coeff[0].show_matrix()
#coeff[1].show_matrix()

print "Verify the orthogonality of orbitals"
C = coeff[0]
S = C.H() * C
S.show_matrix()

print "Orthogonalize orbitals"
Ctilda = utils.orthogonalize_orbitals(C)
S = Ctilda.H() * Ctilda
S.show_matrix()

C = coeff[1]
S = C.H() * C
S.show_matrix()

Ctilda = utils.orthogonalize_orbitals(C)
S = Ctilda.H() * Ctilda
S.show_matrix()




# no SOC, spin-polarized, several kpts
print "================ no SOC, spin-polarized, several kpts ==================="

params = {"nac_method":1, "prefix": "nosoc_spinpol_2x2x2/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
sz = len(e)
for k in xrange(sz):
    print "==== kpt %i ======" % (k)

    print "Energies of the active space orbitals, Ha"
    E = e[k].show_matrix()

    print "Active space orbitals"
    print "uncomment"
    #coeff[k].show_matrix()

    print "Verify the orthogonality of orbitals"
    C = coeff[k]
    S = C.H() * C
    S.show_matrix()


    print "Orthogonalize orbitals"
    Ctilda = utils.orthogonalize_orbitals(C)
    S = Ctilda.H() * Ctilda
    S.show_matrix()




# Method 2
# SOC, non-collinear, one kpt
# Now, all the orbitals are singly-occupied by definition, so the HOMO is = the number of electrons (72).
# But... to keep the parallels with the noSOC cases, we still think in terms of the no-SOC orbitals, that
# is, again, the HOMO will be 36. The transformation to the HOMO=72 case is done inside the function.

print "================ SOC, non-collinear, one kpt ==================="

params = {"nac_method":2, "prefix": "soc_1x1x1/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
E = e[0].show_matrix()

print "Active space orbitals"
print "uncomment"
#coeff[0].show_matrix()
#coeff[1].show_matrix()

print "Verify the orthogonality of orbitals"
C = coeff[0]
S = C.H() * C
S.show_matrix()

print "Orthogonalize orbitals"
Ctilda = utils.orthogonalize_orbitals(C)
S = Ctilda.H() * Ctilda
S.show_matrix()




# SOC, non-collinear, several kpts
# Now, all the orbitals are singly-occupied by definition, so the HOMO is = the number of electrons (72).
# But... to keep the parallels with the noSOC cases, we still think in terms of the no-SOC orbitals, that
# is, again, the HOMO will be 36. The transformation to the HOMO=72 case is done inside the function.
print "================ SOC, non-collinear, several kpts ==================="

params = {"nac_method":2, "prefix": "soc_2x2x2/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
sz = len(e)
for k in xrange(sz):
    print "==== kpt %i ======" % (k)

    print "Energies of the active space orbitals, Ha"
    E = e[k].show_matrix()

    print "Active space orbitals"
    print "uncomment"
    #coeff[k].show_matrix()

    print "Verify the orthogonality of orbitals"
    C = coeff[k]
    S = C.H() * C
    S.show_matrix()

    print "Orthogonalize orbitals"
    Ctilda = utils.orthogonalize_orbitals(C)
    S = Ctilda.H() * Ctilda
    S.show_matrix()


