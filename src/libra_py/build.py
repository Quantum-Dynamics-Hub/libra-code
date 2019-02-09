#*********************************************************************************                     
#* Copyright (C) 2017-2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: build
   :platform: Unix, Windows
   :synopsis: This module implements functions for building molecular structures

.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import LoadPT
import units


def read_xyz(filename, inp_units=1, out_units=0):
    """Read an xyz file (single structure)
    
    Args:
        filename ( string ): the name of the xyz file to read
        inp_units ( int ): defines the coordinates' units in the input file:

            * 0 - Bohr
            * 1 - Angstrom ( default )

        out_units ( int ): defines the output coordinates' units:

            * 0 - Bohr ( default )
            * 1 - Angstrom

    Returns: 
        tuple: (L, coords), where:
 
            * L ( list of N strings ): the labels of N atoms read
            * coords ( list of N VECTORs ): the coordinates of N atoms read [defined by out_units]

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
             
            R = VECTOR(x,y,z) * units.length_converter(inp_units, out_units)

            L.append(tmp[0])
            coords.append(R)

    return L, coords


def make_xyz(L, R, inp_units=0, out_units=1):
    """Convert atomic labels and coordinates to a string formatted according to xyz

    Args:
        L ( list of N strings ): the labels of N atoms in the system
        R ( list of N VECTORs ): the coordinates of N atoms in the system
        inp_units ( int ): defines the units of variables stored in R:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

        out_units ( int ): defines the units of the coordinates written in xyz file:

            * 0 - Bohr 
            * 1 - Angstrom ( default )

    Returns: 
        string: a string representing the xyz file of the system
                

    """

    nat = len(L)
    res = "%i\n" % (nat)
    res = res + "\n"

    for i in xrange(nat):
        res = res + "%s   %12.8f  %12.8f  %12.8f \n" % (L[i], R[i].x/units.Angst, R[i].y/units.Angst, R[i].z/units.Angst)
    return res



def read_xyz_crystal(filename, a,b,c, inp_units=0, out_units=0):
    """Read an xyz file with atomic coordinates given in fractional lattice vector units

    Args:
        filename ( string ): the name of the xyz file to read
        a ( VECTOR ): the unit cell vector in the direction a [units defined by inp_units]
        b ( VECTOR ): the unit cell vector in the direction b [units defined by inp_units]
        c ( VECTOR ): the unit cell vector in the direction c [units defined by inp_units]
        inp_units ( int ): defines the units of a,b,c:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

        out_units ( int ): defines the units of the coordinates in the ourput variable:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

    Returns: 
        tuple: (L, coords), where:
 
            * L ( list of N strings ): the labels of N atoms read
            * coords ( list of N VECTORs ): the coordinates of N atoms read [defined by out_units]
                
    """

    M = MATRIX(3,3)
    M.init(a,b,c)
    M = M.T()

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
             
            R = M * VECTOR(x,y,z)  * units.length_converter(inp_units, out_units)

            L.append(tmp[0])
            coords.append(R)

    return L, coords




def generate_replicas_xyz(tv1, tv2, tv3, rep1, rep2, rep3 , filename, outfile, inp_units=0, out_units=0):
    """

    This function generates a lattice of coordinates and an array of 
    corresponding atomic labels by replicating a given set of labeled
    atoms periodically in up to 3 dimensions (with a variable number of 
    replicas in each direction). This function first reads in the 
    atomic coordinates (supplied in an .xyz format) and then writes 
    the resulting data into another file (.xyz format)

    Args:
        tv1 ( VECTOR ): translation vector in direction 1 [defined by inp_units] 
        tv2 ( VECTOR ): translation vector in direction 2 [defined by inp_units] 
        tv3 ( VECTOR ): translation vector in direction 3 [defined by inp_units] 
        rep1 ( int ): the number of replications along the vector tv1, not counting the original cell
        rep2 ( int ): the number of replications along the vector tv2, not counting the original cell
        rep3 ( int ): the number of replications along the vector tv3, not counting the original cell
        filename ( string ): the name of the file containing input coordinates (original unit cell)
        outfile ( string ): the name of the file where the resulting lattice atoms will be written
        inp_units ( int ): defines the units of variables in filename, tv1, tv2, and tv3:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

        out_units ( int ): defines the units of the coordinates written to outfile

            * 0 - Bohr ( default )
            * 1 - Angstrom 

    Returns:
        None: but a new .xyz file will be printed out

    Example:
        The example below read in a unit cell of Si (8 atoms), replicates in 5 times 
        in x, y, and z directions (so we get a 6 x 6 x 6 supercell). The resulting 
        structure is written to the "cube10.xyz" file

        >>> a = [     5.4307100000000000,    0.0000000000000000,    0.0000000000000000 ]
        >>> b = [     0.0000000000000000,    5.4307100000000000,    0.0000000000000000 ]
        >>> c = [     0.0000000000000000,    0.0000000000000000,    5.4307100000000000 ]
        >>> generate_replicas_xyz(a, b, c, 5, 5, 5 , "Si.xyz", "cube10.xyz")


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



def generate_replicas_xyz2(L, R, tv1, tv2, tv3, Nx, Ny, Nz, inp_units=0, out_units=0):
    """

    This function generates a lattice of coordinates and an array of 
    corresponding atomic labels by replicating a given set of labeled
    atoms periodically in up to 3 dimensions (with a variable number of 
    replicas in each direction)

    Args:
        L ( list of strings ): atomis labels
        R ( list of VECTOR objects ): initial atomic coords [defined by inp_units] 
        tv1 ( VECTOR ): translation vector in direction 1 [defined by inp_units] 
        tv2 ( VECTOR ): translation vector in direction 2 [defined by inp_units] 
        tv3 ( VECTOR ): translation vector in direction 3 [defined by inp_units] 
        Nx ( int ): the number of replications along the vector tv1, not counting the original cell
        Ny ( int ): the number of replications along the vector tv2, not counting the original cell
        Nz ( int ): the number of replications along the vector tv3, not counting the original cell
        inp_units ( int ): defines the units of variables R, tv1, tv2, and tv3:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

        out_units ( int ): defines the units of the coordinates returned:

            * 0 - Bohr ( default )
            * 1 - Angstrom 

    Returns:
        tuple: (lab, newR), where:
 
            * lab ( list of N strings ): the labels of all atoms and their corresponding replicas
            * newR ( list of N VECTORs ): the coordinates of all replicated atoms [defined by out_units]


    """

    nat = len(R)

    lab, newR = [], []
    scl = units.length_converter(inp_units, out_units)

    for i in xrange(nat):

        for nx in range(Nx):
            for ny in range(Ny):
                for nz in range(Nz):
                 
                    # Now apply PBC
                    r = VECTOR(R[i] + nx * tv1 + ny * tv2 + nz * tv3) * scl

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


def crop_sphere_xyz3(L, R, Rcut, pairs, new_L):
    """
    L - list of element names
    R - list of VECTOR objects with the coordinates of the particles, in Bohrs
    Rcut - the radius of the sphere from the center of the QD, in Bohrs
    pairs - list [int, int, VECTOR, VECTOR] - integers describe the indices for the connected atoms
    new_L - is a map containing pairs (string:string), where the key string corresponds to the label of 
          the atom that is inside the cropped sphere, and the value string corresponds to the label of 
          the atom which will turn out to be connected to this atom, but is outside the sphere

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


    # Now add a number of atoms that are ourside the sphere
    added = [] # indices of the added atoms - we need them to avoid a multiple placing
    for it in pairs:
           
        if (it[0] in indx and it[1] in indx):  # both atoms are within the sphere
            pass
        elif (it[0] in indx and it[1] not in indx):  # it[0] is inside, it[1] is outside
            if it[1] not in added:

                new_lab = new_L[ L[it[0]] ]  # new label
                if new_lab == "none":
                    pass
                else:
                    coords.append(R[it[1]])
                    lab.append(new_lab)            
                    added.append(it[1])

        elif (it[1] in indx and it[0] not in indx):  # it[1] is inside, it[0] is outside
            if it[0] not in added:
                new_lab = new_L[ L[it[1]] ]  # new label
                if new_lab == "none":
                    pass
                else:
                    coords.append(R[it[0]])
                    lab.append(new_lab)
                    added.append(it[0])


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




def add_atom_to_system(syst, coords, MaxCoords, Nx,Ny,Nz, a,b,c, shift, elt, mass, scl, max_coord, rnd, inp_units=0):
    """
       
    This function adds an atom and its periordic replicas to a (chemical) System object. 
    The momenta of the atoms are also initialized and are drawn from the normal distribution along
    each of the directions (Cartesian x, y, and z), no center of mass momentum subtraction is made however.
    This is good for constructing lattices.

    Args:
        syst ( System object ): the chemical system to which we are going to add the atoms
        coords ( list of VECTORs ): coordinates of all atoms in a system
        MaxCoords (list of ints): Maximal coordination number of each atom
        Nx ( int ): how many times to replicate the atom along a direction
        Ny ( int ): how many times to replicate the atom along b direction
        Nz ( int ): how many times to replicate the atom along c direction
        a ( VECTOR ): the periodicity translation vector along the a direction [Bohr]
        b ( VECTOR ): the periodicity translation vector along the b direction [Bohr]
        c ( VECTOR ): the periodicity translation vector along the c direction [Bohr]
        shift ( VECTOR ): essentially a coordinate of the atom in the unit cell [Bohr]
        elt ( string ): Element name
        mass ( double ): mass of the atom [a.u. of mass, 1 = mass of an electron]        
        scl ( double ): width of the momentum distrubution [a.u. of momentum]
        max_coord ( int ): maximal coordination number of the added atom allowed
        rnd ( Random ): an instance of the random number generator class
        inp_units ( int ): defines the units of variables stored in shift, a,b, and c

            * 0 - Bohr ( default )
            * 1 - Angstrom 


    Returns:
        None: But the following objects will be updated: 
     
        * syst - the new atoms will be added
        * coords - the new coordinates will be added
        * MaxCoords - the maximal coordination numbers for new atoms will be added

    """

    nat = syst.Number_of_atoms

    # Associate the coordinates with a molecular system
    U = Universe() 
    LoadPT.Load_PT(U, "elements.dat", 0)

    i = 0

    out_units = 0 # internal units in System objects are atomic units
    for nx in xrange(Nx):
        for ny in xrange(Ny):
            for nz in xrange(Nz):
             
                R = VECTOR(nx * a + ny * b + nz * c + shift) * units.length_converter(inp_units, out_units)
                coords.append(R)
                MaxCoords.append(max_coord)

                syst.CREATE_ATOM( Atom(U,  {"Atom_element":elt,"Atom_cm_x":R.x,"Atom_cm_y":R.y,"Atom_cm_z":R.z})  )

                # Masses
                syst.Atoms[nat+i].Atom_RB.set_mass(mass)

                # Momenta
                p = VECTOR(scl.x*rnd.normal(), scl.y*rnd.normal(), scl.z*rnd.normal())
                syst.Atoms[i].Atom_RB.set_momentum(p);

                i = i + 1



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
