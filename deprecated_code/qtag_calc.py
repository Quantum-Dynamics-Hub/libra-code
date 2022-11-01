#*********************************************************************************
#* Copyright (C) 2021-2022 Matthew Dutra, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
"""
..module:: qtag_calc
  :platform: Unix, Windows
  :synopsis: This module contains "ground-level" functions for calculations, mostly output things like energy, norm, etc.

..moduleauthors :: Matthew Dutra, Alexey Akimov
"""

import os, sys
import numpy as np

from liblibra_core import *
from libra_py import data_outs

# CMATRIX qtag_psi(MATRIX& q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff);
def psi(ndof, ntraj_on_surf, qpas, c, x0):
    """Returns the (complex) wavefunction value *wf* at a given point *x0*, calculated using the single-surface basis parameters stored in *qpas* and coefficients *c*.

    Args:
        ndof (integer): Number of degrees of freedom.

        ntraj_on_surf (integer): Number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The ntraj_on_surf-by-1 complex matrix of basis coefficients.

        x0 (MATRIX): The matrix of coordinates [x_1,...,x_ndof] at which the wavefunction value should be calculated.

    Returns:
        wf (complex): Complex value of the wavefunction at x0.
    """

    qvals,pvals,avals,svals=qpas[0],qpas[1],qpas[2],qpas[3]

    wf=0+0j
    for i in range(ntraj_on_surf):
        prod=1.0
        for j in range(ndof):
            q1,p1,a1,s1=qvals.get(j,i),pvals.get(j,i),avals.get(j,i),svals.get(j,i)
            prod*=(a1/np.pi)**0.25*np.exp(-a1/2.0*(x0.get(j)-q1)**2+1j*(p1*(x0.get(j)-q1)+s1))
        wf+=c.get(i)*prod
    return(wf)



def wf_calc_nD(dyn_params, plt_params, prefix):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*,
       based on the input contained in the dict *dyn_params*. The placement is randomly chosen from a
       Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho_cut*,
       and each surface has the same number of trajectories. The corresponding Gaussian-distributed
       momenta are ordered so that the wavepacket spreads as x increases.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`nsteps`]** (int) : the number of simulation steps

          * **dyn_params[`states`]** (int) : the list of states

          * **dyn_params[`grid_dims`]** (list of floats) : the total number of basis functions to be
              placed on each surface. For Gaussian, the list has only one element, specifying the
              number of basis functions per surface. Note that the total number of basis functions
              will then be *nstates*-by-*prod(grid_dims)*

          * **dyn_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

        prefix (str): The name of the directory containing the q, p, a, and s trajectory data.
    """

    #Collect simulation parameters from dyn_params dict...
    nsteps = dyn_params["nsteps"]
    states = dyn_params["states"]
    grid_dims = dyn_params["grid_dims"]
    ndof = dyn_params["ndof"]

    nstates = len(states)
    xmin = plt_params["xmin"]
    xmax = plt_params["xmax"]
    npoints = plt_params["npoints"]

    #This is standard Libra convention for wf files...
    data_type1="wfcr"; data_type2="dens"; data_type3="rep_0"

    #Determine ntraj from dyn_params grid_dims variable...
    ntraj = 1
    for i in range(len(grid_dims)):
        ntraj *= grid_dims[i]

    #Define which timesteps to calculate the wf for...
    lines = plt_params['snaps']

    #Open directory with coeffs, q, p, a, s data...
    qfile = open(prefix+"/q.txt")
    pfile = open(prefix+"/p.txt")
    afile = open(prefix+"/a.txt")
    cfile = open(prefix+"/coeffs.txt")

    #Make the output directory...
    if not os.path.isdir(prefix+"/wfc"):
        os.mkdir(prefix+"/wfc")

    #Create the mesh to calculate the wf on...
    grid_bounds = np.mgrid[tuple(slice(xmin[dof],xmax[dof],complex(0,npoints[dof])) for dof in range(ndof))]

    elems = 1
    for i in npoints:
        elems *= i

    b = grid_bounds.flatten()
    wfpts = []
    for i in range(elems):
        index = i
        coords = []

        while index < len(b):
            coords.append(b[index])
            index += elems
        wfpts.append(coords)

    #Read the trajectory data and compute the wf on the mesh...
    iline=0
    for line in range(nsteps):
        qdata = qfile.readline().strip().split()
        pdata = pfile.readline().strip().split()
        adata = afile.readline().strip().split()
        coeffs = cfile.readline().strip().split()

        if iline in lines:
            outfile = open(prefix+"/wfc/"+data_type1+"_snap_"+str(line)+"_"+data_type2+"_"+data_type3,"w")

            for pt in wfpts:
                for nn in pt:
                    outfile.write(str(nn)+" ")
                idata = 0
                for state in range(nstates):
                    wf = 0+0j
                    for j in range(ntraj):
                        coeff = complex(float(coeffs[2*(j+state*ntraj)]),
                                        float(coeffs[2*(j+state*ntraj)+1]))

                        gaus = 1.0+0.0j
                        for dof in range(ndof):
                            qj = float(qdata[idata])
                            pj = float(pdata[idata])
                            aj = float(adata[idata])
                            idata +=1

                            val = pt[dof]
                            gaus*=(aj/np.pi)**0.25*np.exp(-aj/2.0* \
                                (val-qj)**2+1j*(pj*(val-qj)))

                        wf+=coeff*gaus
                    outfile.write(str(abs(wf)**2)+" ")
                    #outfile.write(str(np.imag(wf))+" ")
                outfile.write("\n")
        iline+=1

