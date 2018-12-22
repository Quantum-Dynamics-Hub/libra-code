#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



def calculate(energy_prefix,energy_suffix,dip_prefix,dip_suffix,isnap,fsnap,opt,scl1,scl2,outfile,HOMO,minE,maxE,dE):
    """
    This fucntions Computes average value of the Hamiltonian matrix elements in energy axis

    \param[in] energy_prefix Prefix of the files containing the energy of the states (diagonal terms)
    \param[in] energy_suffix Similarly to energy_prefix
    \param[in] dip_prefix The prefix of the series of filenames containing transition dipole moment
    \param[in] dip_suffix The suffix of the series of filenames containing transition dipole moment
    \param[in] isnap Index of initial file
    \param[in] fsnap Index of final file
    \param[in] opt Option for averaging
    \param[in] scl1 Scaling factor for energy scale
    \param[in] scl2 Scaling factor for output quantity scale 
    \param[in] outfile The name of the file, where the map will be written
    \param[in] HOMO The index (starting from 1) of the HOMO orbital
    \param[in] inE, maxE, dE Determine the range and the resolution of the absorption spectrum
  
    Returns: A list of grid points for the x-axis and an intensity of unit energy for each
             point on the x-axis   
    """

    filename = energy_prefix + str(isnap) + energy_suffix
    f = open(filename,"r")
    a = f.readline()
    f.close()
    sz = len(a.split())/2
    print "Map size is %d by %d" % (sz,sz)

    HOMO = HOMO - 1 # now it has the meaning of the index rather than id

    #========= Prepare storage ============
    M = []
    i = 0
    while i<sz:
        m = []
        j = 0
        while j<sz:
            m.append(0.0)
            j = j + 1
        M.append(m)
        i = i + 1

    exE = []
    exI = []
    e = minE
    npt = 0
    while e<(maxE+dE):
        exE.append(e)
        exI.append(0.0)
        e = e + dE
        npt = npt + 1


    #========== Read in data ==============
    count = 0.0

    for k in range(isnap,fsnap+1):
        filename = dip_prefix + str(k) + dip_suffix
        e_filename = energy_prefix + str(k) + energy_suffix
        print "Reading files ",filename, e_filename


        if(os.path.exists(filename) and os.path.exists(e_filename)):

            # Read energies of the states
            f1 = open(e_filename,"r")
            A1 = f1.readlines()
            f1.close()

            E = []
            i = 0
            while i<sz:
                tmp1 = A1[i].split()
                y = scl1*float(tmp1[i])
                E.append(y)
                i = i + 1

            # Read transition dipole moments
            f = open(filename,"r")
            A = f.readlines()
            f.close()
            i = 0
            while i<sz:
                tmp = A[i].split()
                j = 0
                while j<sz:
                    x = float(tmp[j])
                    if opt==0:
                        M[i][j] = x
                    elif opt==1:
                        M[i][j] = math.fabs(x)
                    elif opt==2:
                        M[i][j] = x*x

                    j = j + 1
                i = i + 1

            # Now compute transition intensities
            i = HOMO
            while i>0:
                j = HOMO+1
                while j<sz:
                    # ====== Now we are at the point to consider i->j transition
                    de = E[j] - E[i]
                    if(de>minE and de<maxE):
                        indx = int( (de - minE)/dE )
                        exI[indx] = exI[indx] + M[i][j]
                        count = count + 1.0

                    j = j + 1
                i = i-1


    #========== Average and print out results ===============
    out = open(outfile,"w")

    exI[0] = 0.0   # Do not print dipole moments
    i = 0
    while i<npt:
        if(opt==0 or opt==1):
            exI[i] = scl2*exI[i] / count
        elif opt==2:
            exI[i] = scl2 * math.sqrt(exI[i]) / count

        out.write("%f  %f \n" % (exE[i], exI[i]) )
        line = str(exE[i]) + "  "+str(exI[i]) + "\n"
        out.write(line)
        i = i + 1

    out.close()
    return exE, exI



def ham_map(prefix, isnap, fsnap, suffix, opt, scl, outfile):
    """
    This fucntions Computes the average value of the Hamiltonian matrix elements in energy axis:

    \param[in] prefix The prefix of the series of filenames containing required matrix elements 
    \param[in] isnap Index of initial file
    \param[in] fsnap Index of final file
    \param[in] suffix Similarly to prefix
    \param[in] opt Option for averaging
    \param[in] scl Scaling factor. e.g. to convert units
    \param[in] outfile The name of the file, where the map will be written
    """

    filename = prefix + str(isnap) + suffix
    f = open(filename,"r")
    a = f.readline()
    f.close()
    sz = len(a.split())
    print "Map size is %d by %d" % (sz,sz)

    #========= Prepare storage ============
    M = []
    i = 0
    while i<sz:
        m = []
        j = 0
        while j<sz:
            m.append(0.0)
            j = j + 1
        M.append(m)
        i = i + 1

    #========== Read in data ==============
    count = 0.0

    for k in range(isnap,fsnap+1):
        filename = prefix + str(k) + suffix
        print "Reading file ",filename

        if(os.path.exists(filename)):

            f = open(filename,"r")
            A = f.readlines()
            f.close()

            i = 0
            while i<sz:
                tmp = A[i].split()
                j = 0
                while j<sz:
                    x = float(tmp[j])
                    if opt==0:
                        M[i][j] = M[i][j] + x
                    elif opt==1:
                        M[i][j] = M[i][j] + math.fabs(x)
                    elif opt==2:
                        M[i][j] = M[i][j] + x*x

                    j = j + 1
                i = i + 1

            count = count + 1.0

    #========== Average and print results ===============
    out = open(outfile,"w")

    i = 0
    while i<sz:
        j = 0
        while j<sz:
            # ====== Average ======
            if(opt==0 or opt==1):
                M[i][j] = scl * M[i][j] / count
            elif opt==2:
                M[i][j] = scl * math.sqrt(M[i][j] / count )

            # ======= Print =======            
            if(i==j):
                M[i][j] = 0.0  # because diagonal and off-diagonal elements are of different orders of magnitude
            out.write("%d  %d  %f\n" % (i, j, M[i][j]))
            j = j + 1
        out.write("\n")
        i = i + 1


    return 1



