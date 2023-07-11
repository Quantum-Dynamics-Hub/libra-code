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


#===================== Reading info ==================================
# Method 0 
# no SOC, non-polarized, one kpt
# To get HOMO, divide the number of electrons (72) by 2 => 36. Everything else is relative to it
print "================ no SOC, non-polarized, one kpt ==================="

params = {"nac_method":0, "prefix": "nosoc_nospinpol_1x1x1/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
Edia = e[0]; e[0].show_matrix()

print "Verify the orthogonality of orbitals"
Cdia = coeff[0]
S = Cdia.H() * Cdia
S.show_matrix()



# Method 1
# no SOC, spin-polarized, one kpt
# To get HOMO, divide the number of electrons (72) by 2 => 36. Everything else is relative to it
print "================ no SOC, spin-polarized, one kpt ==================="

params = {"nac_method":1, "prefix": "nosoc_spinpol_1x1x1/x0.export", "minband":36, "maxband":38}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
Edia_a = e[0]; e[0].show_matrix()
Edia_b = e[1]; e[1].show_matrix()

print "Verify the orthogonality of orbitals"
Cdia_a = coeff[0]
S = Cdia_a.H() * Cdia_a
S.show_matrix()

Cdia_b = coeff[1]
S = Cdia_b.H() * Cdia_b
S.show_matrix()



# Method 2
# SOC, non-collinear, one kpt
# Now, all the orbitals are singly-occupied by definition, so the HOMO is = the number of electrons (72).
# But... to keep the parallels with the noSOC cases, we still think in terms of the no-SOC orbitals, that
# is, again, the HOMO will be 36. The transformation to the HOMO=72 case is done inside the function.

print "================ SOC, non-collinear, one kpt ==================="

params = {"nac_method":2, "prefix": "soc_1x1x1/x0.export", "minband":1, "maxband":50}
info, e, coeff, grid = trajectory.read_all(params)

print info, e, coeff, grid
print "Energies of the active space orbitals, Ha"
Eadi = e[0]; e[0].show_matrix()

print "Verify the orthogonality of orbitals"
Cadi = coeff[0]
S = Cadi.H() * Cadi
S.show_matrix()


print "================ PROJECTIONS ==================="
print "Projections of SOC wfc. onto non-spin-polarized ones and vice-versa"

L = Cdia.H() * Cadi
print "<Psi_dia|Psi_adi>"
L.show_matrix()

Lt = Cadi.H() * Cdia
print "<Psi_adi|Psi_dia>"
Lt.show_matrix()


print "Projections of SOC wfc. onto spin-polarized ones and vice-versa"

L = Cdia_a.H() * Cadi
print "<Psi_dia(alp)|Psi_adi>"
L.show_matrix()

L = Cdia_b.H() * Cadi
print "<Psi_dia(bet)|Psi_adi>"
L.show_matrix()


Lt = Cadi.H() * Cdia_a
print "<Psi_adi|Psi_dia(alp)>"
Lt.show_matrix()

Lt = Cadi.H() * Cdia_b
print "<Psi_adi|Psi_dia(bet)>"
Lt.show_matrix()


