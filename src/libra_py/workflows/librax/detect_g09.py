#*********************************************************************************
#* Copyright (C) 2017 Kosuke Sato, Alexey V. Akimov
#* Olga S. Bokareva
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file detect_g09.py
# This module implements a function that detect columns showing
# gradients, molecular energies, molecular orbitals, and atomic basis information
# written in the Gaussian09 output file; also implements a function that shows 
# the detected columns for debugging.

import os, re
import sys
import math


def detect_columns(inp_lines):
    ##
    # Finds the keywords and their patterns and extracts the descriptors info
    # \param[in] inp_lines The list of lines containing G09 output file to be unpacked
    # 
    # info - The returned dictionary of descriptors of the given input lines.
    #
    # Used in: detect_g09.py/detect

    info = {}

    i = -1
    
    for line in inp_lines:
        i += 1
        spline = line.split()

        # the number of electrons
        # in G09, alpha and beta electrons are listed separately, so we add up them
        if len(spline) == 6 and spline[2] == "electrons":
	    print "done with electrons"
            info["lele"] = i
            info["Nele"] = int(spline[0])+int(spline[3])

        # the number of atoms is not given in G09-out file. That is why is it counted in separate python loop
        #if len(spline) == 6 and spline[0] == "TOTAL" and spline[3]== "ATOMS":
        #    info["latoms"] = i
        #    info["Natoms"] = int(spline[5])

        # the number of occupied orbitals (alpha and beta) is given in one line
        if len(spline) == 6 and spline[1] == "alpha" and spline[4] == "beta":
    	    print "done with alpha and beta"
            info["locc_alp"] = i
            info["Nocc_alp"] = int(spline[0])
            info["locc_bet"] = i
            info["Nocc_bet"] = int(spline[3])

        # the number of cartesian gaussian basis functions
        if len(spline) == 10 and "cartesian basis functions" in line:
            print "done with cartesian basis func"
            info["lgbf"] = i
            info["Ngbf"] = int(spline[6])

        # the atomic basis sets (normal basis sets!), no semiempirical or ECP so far
        if "AO basis set in the form of general basis input" in line:
    	    info["ab_start"] = i + 1
    	    print "done with ab_start"
        if "AO basis set (Overlap normalization):" in line:
            print "done with ab_end"
            info["ab_end"] = i - 3
        
        # eigenvectors
        if "Molecular Orbital Coefficients" in line:
    	    print "done with MO"
            info["mo_start"] = i + 1
            info["mo_end"] = i + (info["Ngbf"] + 3) * int(math.ceil(info["Ngbf"]/5.0))

        # the coordinates of the atoms (in Angstroms!!!!) will be given simulteneously with atom counting section
        #if len(spline) == 4 and spline[2] == "COORDINATES" and spline[3] == "(BOHR)":
        #    info["coor_start"] = i + 2

        # the gradients(in Hartree/Bohr)

        if len(spline) == 4 and spline[2] == "Forces" and spline[3] == "(Hartrees/Bohr)":
            print "done with forces"
            info["grad_start"] = i + 3

        # total energy
        # SCF Done:  E(RPM3) = -0.846886626509E-01 A.U. after    9 cycles
#	if "SCF Done:" in line:
        if len(spline) > 5 and  spline[1] == "Done:": # and spline[8] == "cycles":
            print "done with energy"
    	    info["ltot_ene"] = i
            info["tot_ene"] = float(spline[4])

    # separate block for counting atoms in molecule and reading the lines with coordinates
    found = False
    block = []
    i=-1
    for line in inp_lines:
        i += 1
        if found:
	    line=line.strip()
	    block.append(line)

    	    if "Distance matrix" in line: 
    		break
    		info["latoms"] = i
	else:
    	    if "Input orientation" in line:
		info["coor_start"] = i+5
        	found = True
    		block = []

    Natoms = 0
    for i in block:
	elements = i.split()
	if len(elements) == 6 and re.search('[a-zA-Z]+',i) == None:
	    Natoms += 1
    info["Natoms"] = Natoms

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
#    print "according to the %i th column," % (info["latoms"]+1)
#    print inp_lines[info["latoms"]]
    print "the number of atoms was counted"
    print "Natoms = %i" % info["Natoms"]
    print "*******************************************"
    print 
    print "******************************************"
    print "according to the %i th column," % (info["lele"]+1)
    print inp_lines[info["lele"]].strip()
    print "Nele = %i" % info["Nele"]
    print "*******************************************"
    print 
    print "******************************************"
    print "according to the %i th column," % (info["locc_alp"]+1)
    print inp_lines[info["locc_alp"]].strip()
    print "Nocc_alp = %i" % info["Nocc_alp"]
    print "according to the %i th column," % (info["locc_bet"]+1)
    print inp_lines[info["locc_bet"]].strip()
    print "Nocc_bet = %i" % info["Nocc_bet"]
    print "*******************************************"
    print
    print "******************************************"
    print "according to the %i th column," % (info["lgbf"]+1)
    print inp_lines[info["lgbf"]].strip()
    print "Ngbf = %i" % info["Ngbf"]
    print "*******************************************"
    print

    #if flag_ao == 1:
    print "******************************************"
    print "ATOMIC BASIS SET is within %i - %i th lines." % (info["ab_start"]+1, info["ab_end"]+1)
    for i in range(info["ab_start"],info["ab_end"]+1):
        print inp_lines[i].strip()
    print "******************************************"
    print

    print "******************************************"
    print "MOLECULAR ORBITALS is within %i - %i th lines." % (info["mo_start"]+1,info["mo_end"]+1)
    for i in range(info["mo_start"],info["mo_end"]+1):
        print inp_lines[i].strip()
    print "******************************************"
    print
    print "******************************************"
    print "COORDINATES OF ATOMS (in Bohr) is within %i - %i th lines." %(info["coor_start"]+1,info["coor_end"]+1)
    for i in range(info["coor_start"],info["coor_end"]+1):
        print inp_lines[i].strip()
    print "******************************************"
    print
    print "******************************************"
    print "GRADIENT (in Hartree/Bohr) is within %i - %i th lines." %(info["grad_start"]+1,info["grad_end"]+1)
    for i in range(info["grad_start"],info["grad_end"]+1):
        print inp_lines[i].strip()
    print "******************************************"
    print
    print "******************************************"
    print "according to the %i th column," % (info["ltot_ene"]+1)
    print inp_lines[info["ltot_ene"]].strip()
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


