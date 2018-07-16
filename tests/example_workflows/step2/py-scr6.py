#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/
from PYXAID2 import *
import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

#cwd = os.getcwd()
#print "Current working directory", cwd
sys.path.insert(1,"/projects/academic/alexeyak/brendan/libra/libra-code/workflows/pyxaid2")

import excitation_spectrum


# Example of usage
opt = 1
scl1 = 13.60569253 # Ry to eV
scl2 = 1000.0 # some scaling to get plottable data
tmin = 0
tmax = 1500

HOMO = 11
minE = 0.0
maxE = 10.0
dE = 0.05


# Energy, couplings and H' in space of orbital indexes
#excitation_spectrum.ham_map("res/hvib_adi_",   tmin,tmax,"_re" ,opt,scl1,"spectr/ave_Ham_re.dat")
#excitation_spectrum.ham_map("res/hvib_adi_",   tmin,tmax,"_im" ,opt,scl1,"spectr/ave_Ham_im.dat")
excitation_spectrum.ham_map("res/hvib_dia_",   tmin,tmax,"_re" ,opt,scl1,"spectr/ave_Ham_re.dat")
excitation_spectrum.ham_map("res/hvib_dia_",   tmin,tmax,"_im" ,opt,scl1,"spectr/ave_Ham_im.dat")


# Same, but in space of orbital energies
#utils.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,scl1,"spectr/1ave_Ham_re.dat")
#utils.ham_map1("res/0_Ham_","_re","res/0_Ham_", tmin,tmax,"_im" ,opt,scl1,scl1,"spectr/1ave_Ham_im.dat")



# Absorption Spectrum
#[Ex, Ix] = excitation_spectrum.calculate("res/E_adi_ks_","_re","res/0_Hprime_","x_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectrx.dat",HOMO,minE,maxE,dE)
#[Ey, Iy] = excitation_spectrum.calculate("res/E_adi_ks_","_re","res/0_Hprime_","y_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectry.dat",HOMO,minE,maxE,dE)
#[Ez, Iz] = excitation_spectrum.calculate("res/E_adi_ks_","_re","res/0_Hprime_","z_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectrz.dat",HOMO,minE,maxE,dE)
"""
[Ex, Ix] = excitation_spectrum.calculate("res/E_dia_ks_","_re","res/0_Hprime_","x_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectrx.dat",HOMO,minE,maxE,dE)
[Ey, Iy] = excitation_spectrum.calculate("res/E_dia_ks_","_re","res/0_Hprime_","y_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectry.dat",HOMO,minE,maxE,dE)
[Ez, Iz] = excitation_spectrum.calculate("res/E_dia_ks_","_re","res/0_Hprime_","z_im",tmin,tmax,opt,scl1,scl2,"spectr/ab_spectrz.dat",HOMO,minE,maxE,dE)


#===================== User scripting goes here ===================

# Convert the 'naked' data in to the format for convolution
XP, YP, ZP = [], [], []
i = 0
while i<len(Ix):
    XP.append([Ex[i], Ix[i]])
    YP.append([Ey[i], Iy[i]])
    ZP.append([Ez[i], Iz[i]])
    i = i + 1

# Now do the convolution
dx0 = dE  # original grid spacing
dx = dE  # this is new grid spacing
var = 2.0*dE  # new variance

#print XP

TMx,Rx = pdos.convolve(XP,dx0,dx,var)
TMy,Ry = pdos.convolve(YP,dx0,dx,var)
TMz,Rz = pdos.convolve(ZP,dx0,dx,var)

#print R
pdos.printout("spectr/ab_spectrx_conv.dat",TMx,2,Rx)
pdos.printout("spectr/ab_spectry_conv.dat",TMy,2,Ry)
pdos.printout("spectr/ab_spectrz_conv.dat",TMz,2,Rz)
"""
