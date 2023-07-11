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
params = {"nac_method":0, "wd0":"nosoc_nospinpol_1x1x1/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info

# no SOC, non-polarized, several kpts
params = {"nac_method":0, "wd0":"nosoc_nospinpol_2x2x2/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info


# Method 1
# no SOC, polarized, one kpt
params = {"nac_method":1, "wd0":"nosoc_spinpol_1x1x1/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info

# no SOC, polarized, several kpts
params = {"nac_method":1, "wd0":"nosoc_spinpol_2x2x2/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info


# Method 2
# SOC, non-collinear, one kpt
params = {"nac_method":2, "wd1":"soc_1x1x1/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info

# SOC, non-collinear, several kpts
params = {"nac_method":2, "wd1":"soc_2x2x2/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info



# Method 3
# no SOC, polarized, one kpt
# SOC, non-collinear, one kpt
params = {"nac_method":3, "wd0":"nosoc_spinpol_1x1x1/x0.export", "wd1":"soc_1x1x1/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info

# no SOC, polarized, several kpts
# SOC, non-collinear, several kpts
params = {"nac_method":3, "wd0":"nosoc_spinpol_2x2x2/x0.export", "wd1":"soc_2x2x2/x0.export"}
info0, all_e_dum0, info1, all_e_dum1 = trajectory.read_info(params)
print info0  # is a dictionary with useful volumetric info
print info1  # is a dictionary with useful volumetric info



