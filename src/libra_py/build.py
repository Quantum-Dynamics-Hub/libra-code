#*********************************************************************************                     
#* Copyright (C) 2017 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file build.py 
# This module implements functions for building molecular structures
#


import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import LoadPT


Angst_to_Bohr = 1.889725989  

def read_xyz(filename):
    """
    And convert to Bohrs
    """

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    tmp = A[0].split()
    Natoms = int(float(tmp[0]))

    L = [] 
    coords = []

    for a in A[2:]:
        tmp = a.split()
        if len(tmp)==4:
            x = float(tmp[1]) 
            y = float(tmp[2])
            z = float(tmp[3])
             
            # Now apply PBC
            R = VECTOR(x,y,z) * Angst_to_Bohr

            L.append(tmp[0])
            coords.append(R)

    return L, coords



def generate_replicas_xyz2(L, R, tv1, tv2, tv3, Nx, Ny, Nz):
    """
    L - initial labels (list of strings)
    R - initial coords (list of VECTOR objects) in Bohrs
    tv1, tv2, tv3 - translation vectors (VECTOR objects), in Bohrs
    Nx, Ny, Nz - the number of replications along each vector
    """

    nat = len(R)

    lab, newR = [], []

    for i in xrange(nat):

        for nx in range(Nx):
            for ny in range(Ny):
                for nz in range(Nz):
                 
                    # Now apply PBC
                    r = VECTOR(R[i] + nx * tv1 + ny * tv2 + nz * tv3)

                    lab.append(L[i])
                    newR.append(r)

    return lab, newR




def crop_sphere_xyz2(L, R, Rcut):
    """
    L - list of element names
    R - list of VECTOR objects with the coordinates of the particles, in Bohrs
    Rcut - the radius of the sphere from the center of the QD, in Bohrs

    Output:
    new lists of L and R after cropping

    """


    # Compute COM
    Rcom = VECTOR(0.0, 0.0, 0.0)

    for r in R:
        Rcom = Rcom + r

    # Geometric center
    Natoms = len(R)
    Rcom = Rcom / Natoms

 
    # Compute remaining number of atoms
    Nat = 0
    indx = []
    for i in range(0,Natoms):
        r = (R[i] - Rcom).length()
        if(r<=Rcut):
            Nat = Nat + 1
            indx.append(i)
            

    # Prepare the coordinates and labels of the remaining atoms
    coords = []
    lab = []     
    for i in indx:
        coords.append(R[i])
        lab.append(L[i])

    return lab, coords




def add_atoms_to_system(syst, L, R, scl, mass, univ_file):
    """
    syst - the System object, initially may be empty
    L - list of element names
    R - list of VECTOR objects with the coordinates of the particles - in Bohrs
    scl - a VECTOR - scaling coefficients for sampling the initial momenta
    mass - list of atomic masses (floats)
    univ_file - name of the file that contains the Universe properties
    """

    rnd = Random()

    nat = syst.Number_of_atoms

    # Associate the coordinates with a molecular system
    U = Universe() 
    LoadPT.Load_PT(U, univ_file, 0)

    sz = len(R)
    for i in xrange(sz):

        syst.CREATE_ATOM( Atom(U,  {"Atom_element":L[i],"Atom_cm_x":R[i].x,"Atom_cm_y":R[i].y,"Atom_cm_z":R[i].z})  )

        # Masses
        syst.Atoms[nat+i].Atom_RB.set_mass(mass[i])

        # Momenta
        p = VECTOR(scl.x*rnd.normal(), scl.y*rnd.normal(), scl.z*rnd.normal())
        syst.Atoms[nat+i].Atom_RB.set_momentum(p);




def generate_replicas_xyz(tv1,tv2,tv3, rep1, rep2, rep3 , filename, outfile):
    """
    tv1, tv2, tv3 - translation vectors (lists of 3 floats)
    rep1, rep2, rep3 - the number of replications along each vector
    filename - the name of the file (xyz format) containing the unit cell structure
    outfile - the name of the output file
    """


    transl = []
    for x in range(0,rep1+1):
        for y in range(0,rep2+1):
            for z in range(0,rep3+1):
                transl.append([x,y,z])


    f = open(filename,"r")
    A = f.readlines()
    f.close()

    tmp = A[0].split()
    Natoms = int(float(tmp[0]))
    f = open(outfile, "w")
    f.write("%i \n" % int(len(transl) * Natoms))
    f.write(A[1][:-1] + "\n")

    L = [] 
    coords = []

    for T in transl:
        for a in A[2:]:
            tmp = a.split()
            if len(tmp)==4:
                x = float(tmp[1]) 
                y = float(tmp[2])
                z = float(tmp[3])
                 
                # Now apply PBC
                dx = T[0]*tv1[0] + T[1]*tv2[0] + T[2]*tv3[0]
                dy = T[0]*tv1[1] + T[1]*tv2[1] + T[2]*tv3[1]
                dz = T[0]*tv1[2] + T[1]*tv2[2] + T[2]*tv3[2]

                f.write("%s  %8.5f  %8.5f  %8.5f \n" % (tmp[0], x+dx, y+dy, z+dz ) )

                L.append(tmp[0])
                coords.append(VECTOR(x+dx, y+dy, z+dz))
    f.close()

    return L, coords


def crop_sphere_xyz(infile, outfile, Rcut):
    """
    infile - the name of the input file with the input structure (xyz format)
    outfile - the name of the output file where the resulting structure be printed out (xyz format)
    """

    f = open(infile,"r")
    A = f.readlines()
    f.close()

    tmp = A[0].split()
    Natoms = int(float(tmp[0]))
    print A[1][:-1]

    # Read coordinates and compute COM
    L = []
    R = []
    Rcom = [0.0,0.0,0.0]

    for a in A[2:]:
        tmp = a.split()
        if len(tmp)==4:
            x = float(tmp[1]) 
            y = float(tmp[2])
            z = float(tmp[3])

            Rcom[0] = Rcom[0] + x
            Rcom[1] = Rcom[1] + y
            Rcom[2] = Rcom[2] + z

            L.append(tmp[0])
            R.append([x,y,z])

    # Geometric center
    Rcom[0] = Rcom[0] / Natoms
    Rcom[1] = Rcom[1] / Natoms
    Rcom[2] = Rcom[2] / Natoms

 
    # Compute remaining number of atoms
    Nat = 0
    indx = []
    for i in range(0,Natoms):
        r = math.sqrt((R[i][0] - Rcom[0])**2 + (R[i][1] - Rcom[1])**2 + (R[i][2] - Rcom[2])**2)
        if(r<=Rcut):
            Nat = Nat + 1
            indx.append(i)
            

    # Print resulting coordinates
    f1 = open(outfile,"w")
    f1.write("%5i\n" % Nat)
    f1.write(A[1])


    coords = []
    lab = []     
    for i in indx:
        coords.append(VECTOR(R[i][0], R[i][1], R[i][2]))
        lab.append(L[i])
        f1.write("%s  %12.6f  %12.6f  %12.6f \n" % (L[i], R[i][0], R[i][1], R[i][2]) )

    f1.close()

    return lab, coords




if __name__=='__main__':
    a = [     5.4307100000000000,    0.0000000000000000,    0.0000000000000000 ]
    b = [     0.0000000000000000,    5.4307100000000000,    0.0000000000000000 ]
    c = [     0.0000000000000000,    0.0000000000000000,    5.4307100000000000 ]
    generate_replicas_xyz(a, b, c, 5, 5, 5 , "Si.xyz", "cube10.xyz")
    crop_sphere_xyz("cube10.xyz", "qd_5.xyz", 5.0)
