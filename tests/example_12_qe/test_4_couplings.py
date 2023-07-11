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

nmin = 30
params = {"nac_method":0, "prefix": "nosoc_nospinpol_1x1x1/x0.export", "minband":nmin, "maxband":40 }
info, e, coeff, grid = trajectory.read_all(params)

Edia = e[0]
#Cdia = coeff[0]
Cdia = utils.orthogonalize_orbitals(coeff[0])


# Method 2
# SOC, non-collinear, one kpt
# Now, all the orbitals are singly-occupied by definition, so the HOMO is = the number of electrons (72).
# But... to keep the parallels with the noSOC cases, we still think in terms of the no-SOC orbitals, that
# is, again, the HOMO will be 36. The transformation to the HOMO=72 case is done inside the function.

print "================ SOC, non-collinear, one kpt ==================="

params = {"nac_method":2, "prefix": "soc_1x1x1/x0.export", "minband":30, "maxband":40 }
info, e, coeff, grid = trajectory.read_all(params)

Eadi = e[0]; 
#Cadi = coeff[0]
Cadi = utils.orthogonalize_orbitals(coeff[0])


print "================ SOC in the diabatic basis ==================="
# H_dia * Cdia = Cdia * E_dia
# H_dia(no soc) = Cdia * E_dia * Cdia.H()

Hdia =  Cdia * Edia * Cdia.H()
Hadi =  Cadi * Eadi * Cadi.H()

Hsoc = Hadi - Hdia  # perturbation in the PW basis
Esoc = Cdia.H() * Hsoc * Cdia  # perturbation in the spin-diabatic basis

print "SO couplings"
esoc = CMATRIX(3,3)
templ = Py2Cpp_int([36-nmin, 36-nmin+1, 36-nmin+2])
pop_submatrix(Esoc, esoc, templ, templ )

esoc.show_matrix()






