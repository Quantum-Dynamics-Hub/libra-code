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

## \file detect_gms.py
# This module implements a function that detect columns showing
# gradients, molecular energies, molecular orbitals, and atomic basis information
# written in the GAMESS output file; also implements a function that shows 
# the detected columns for debugging.

import os
import sys
import math


def detect_columns(inp_lines):
    ##
    # Finds the keywords and their patterns and extracts the descriptors info
    # \param[in] inp_lines The list of lines containing GAMESS output file to be unpacked
    # 
    # info - The returned dictionary of descriptors of the given input lines.
    #
    # Used in: detect_gms.py/detect

    info = {}

    i = -1
    for line in inp_lines:
        i += 1
        spline = line.split()

        # the number of electrons
        if len(spline) == 5 and spline[2] == "ELECTRONS":
            info["lele"] = i
            info["Nele"] = int(spline[4])

        # the number of atoms
        if len(spline) == 6 and spline[0] == "TOTAL" and spline[3]== "ATOMS":
            info["latoms"] = i
            info["Natoms"] = int(spline[5])

        # the number of occupied orbitals (alpha and beta)
        if len(spline) == 7 and spline[4] == "(ALPHA)":
            info["locc_alp"] = i
            info["Nocc_alp"] = int(spline[6])
        if len(spline) == 8 and spline[4] == "(BETA":
            info["locc_bet"] = i
            info["Nocc_bet"] = int(spline[7])

        # the number of cartesian gaussian basis functions
        if len(spline) == 8 and spline[5] == "FUNCTIONS":
            info["lgbf"] = i
            info["Ngbf"] = int(spline[7])

        # the atomic basis sets

        if len(spline) == 3 and spline[1] == "BASIS" and spline[2] == "SET":
            info["ab_start"] = i + 7
        if len(spline) == 8 and spline[5] == "SHELLS":
            info["ab_end"] = i - 2
        
        # eigenvectors
        if len(spline) > 0 and spline[0] == "EIGENVECTORS":
            info["mo_start"] = i + 3
            
            info["mo_end"] = i + 1 + (info["Ngbf"] + 4) * int(math.ceil(info["Ngbf"]/5.0))

        # the coordinates of the atoms (in Bohr)
        if len(spline) == 4 and spline[2] == "COORDINATES" and spline[3] == "(BOHR)":
            info["coor_start"] = i + 2

        # the gradients(in Hartree/Bohr)

        if len(spline) == 4 and spline[0] == "GRADIENT" and spline[3] == "ENERGY":
            info["grad_start"] = i + 4

        # total energy

        if len(spline) == 8 and spline[0] == "FINAL" and spline[2] == "ENERGY":
            info["ltot_ene"] = i
            info["tot_ene"] = float(spline[4])

    # end index of Coordinates and Gradients

    info["coor_end"] = info["Natoms"] + info["coor_start"] - 1
    info["grad_end"] = info["Natoms"] + info["grad_start"] - 1 
    
    return info


def show_outputs(inp_lines,info):
    ## Find the keywords below
    # \param[in] inp_lines The list of lines containing GAMESS output file to be unpacked
    # \param[in] info The dictionary of descriptors of the given input lines.
    #
    # This function shows the positions of the data elements in the analyzed file and 
    # some other auxiliary information extracted from the file
    #
    # Used in: detect_gms.py/detect

    print "******************************************"
    print "according to the %i th column," % (info["latoms"]+1)
    print inp_lines[info["latoms"]]
    print "Natoms = %i" % info["Natoms"]
    print "*******************************************"
    print 
    print "******************************************"
    print "according to the %i th column," % (info["lele"]+1)
    print inp_lines[info["lele"]]
    print "Nele = %i" % info["Nele"]
    print "*******************************************"
    print 
    print "******************************************"
    print "according to the %i th column," % (info["locc_alp"]+1)
    print inp_lines[info["locc_alp"]]
    print "Nocc_alp = %i" % info["Nocc_alp"]
    print "according to the %i th column," % (info["locc_bet"]+1)
    print inp_lines[info["locc_bet"]]
    print "Nocc_bet = %i" % info["Nocc_bet"]
    print "*******************************************"
    print
    print "******************************************"
    print "according to the %i th column," % (info["lgbf"]+1)
    print inp_lines[info["lgbf"]]
    print "Ngbf = %i" % info["Ngbf"]
    print "*******************************************"
    print

    #if flag_ao == 1:
    print "******************************************"
    print "ATOMIC BASIS SET is within %i - %i th lines." % (info["ab_start"]+1, info["ab_end"]+1)
    for i in range(info["ab_start"],info["ab_end"]+1):
        print inp_lines[i]
    print "******************************************"
    print

    print "******************************************"
    print "MOLECULAR ORBITALS is within %i - %i th lines." % (info["mo_start"]+1,info["mo_end"]+1)
    for i in range(info["mo_start"],info["mo_end"]+1):
        print inp_lines[i]
    print "******************************************"
    print
    print "******************************************"
    print "COORDINATES OF ATOMS (in Bohr) is within %i - %i th lines." %(info["coor_start"]+1,info["coor_end"]+1)
    for i in range(info["coor_start"],info["coor_end"]+1):
        print inp_lines[i]
    print "******************************************"
    print
    print "******************************************"
    print "GRADIENT (in Hartree/Bohr) is within %i - %i th lines." %(info["grad_start"]+1,info["grad_end"]+1)
    for i in range(info["grad_start"],info["grad_end"]+1):
        print inp_lines[i]
    print "******************************************"
    print
    print "******************************************"
    print "according to the %i th column," % (info["ltot_ene"]+1)
    print inp_lines[info["ltot_ene"]]
    print "total energy = %f" % info["tot_ene"]
    print "******************************************"
    print


def detect(inp_lines,flag):
    ## 
    # This function detects the positions of the valuable data in a file represented as a  
    # list of lines. It does not return the data itself, only the descriptors of where to
    # get the info about: atomic basis sets, molecular energies , molecular orbitals,
    # atomic gradients, coordinates, and labels.
    # \param[in] inp_lines The list of lines containing the file to be unpacked
    # \param[in] flag Debug info printing: 1 - print, otherwise - don't 
    # info - is the returned dictionary of descriptors of the given input lines
    #
    # Used in: extract_gms.py/gms_extract

    info = detect_columns(inp_lines)
    
    if flag == 1:
        show_outputs(inp_lines,info)

        print "*********************************************"
        print "detect program ends"
        print "*********************************************\n"

    return info


