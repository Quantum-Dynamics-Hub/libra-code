#*********************************************************************************
#* Copyright (C) 2018-2019 Brendan A. Smith, Alexey V. Akimov
#* Copyright (C) 2016-2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: QE_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for dealing with the outputs from QE (Quantum Espresso) package

.. moduleauthor:: 
       Alexey V. Akimov 
       Brendan A. Smith
  
"""


import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import units
import regexlib as rgl

def cryst2cart(a1,a2,a3,r):
    """Crystal to Cartesian coordinate conversion 

    An auxiliary function that convert vector <r> in crystal (alat) 
    coordinates with the cell defined by vectors a1, a2, a3, to the
    Cartesian coordinate xyz

    Args:
        a1 ( [double, double, double] ): one of the unit cell/periodicty vectors
        a2 ( [double, double, double] ): one of the unit cell/periodicty vectors
        a3 ( [double, double, double] ): one of the unit cell/periodicty vectors
        r ( [double, double, double] ): Crystal (fractional) coordinates of the atom


    Returns:
        [double, double, double]: the Cartesian coordinates of an atom

    """
    x = [0.0,0.0,0.0]
    for i in [0,1,2]:
        x[i] = a1[i]*r[0] + a2[i]*r[1] + a3[i]*r[2]

    return x



def read_qe_index(filename, orb_list, verbose=0):
    """

    This functions reads an ASCII/XML format file containing wavefunction
    and returns the coefficients of the plane waves that constitute the
    wavefunction

    Args:
        filename ( string ): This is the name of the file we will be reading to construct a wavefunction
        orb_list ( list of ints ): The indices of the orbitals which we want to consider. Orbitals indexing 
            at 1, not 0
        verbose ( int ): The flag controlling the amout of extra output:
        
            * 0 - no extra output (default)
            * 1 - print extra stuff
  
    Returns: 
        CMATRIX(norb, norb): a diagonal matrix (complex-valued) of orbital energies - only for those that are of interest
            Here, norb = len(orb_list) and the matrix will only contain numbers that correspond to orbitals defined
            in the orb_list argument [units: Ha]
    """

    Ry2Ha = 0.5 # conversion factor
   
    ctx = Context(filename)  #("x.export/index.xml")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children("Root/Eigenvalues")  #("Kpoint.1")

    info = {}
    info["nspin"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nspin","-1.0")))  # 1 - non-polarized, 2 - polarized, 4 - non-collinear
    info["nk"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nk","-1.0")))  # the number of k-points
    info["nbnd"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nbnd","-1.0")))  # the number of orbitals in each k-point
    info["efermi"] = Ry2Ha * float(ctx.get("Eigenvalues/<xmlattr>/efermi","-1.0"))  # Fermi energy

    # Usially in atomic units!
    info["alat"] = float(ctx.get("Cell/Data/<xmlattr>/alat","-1.0"))    # lattice constant (a)
    info["omega"] = float(ctx.get("Cell/Data/<xmlattr>/omega","-1.0"))  # unit cell volume
    info["tpiba"] = float(ctx.get("Cell/Data/<xmlattr>/tpiba","-1.0"))  # 2 pi / a
    info["tpiba2"] = float(ctx.get("Cell/Data/<xmlattr>/tpiba2","-1.0"))  # ( 2 pi / a )**2

    # Direct lattice vectors
    tmp = ctx.get("Cell/a1/<xmlattr>/xyz","-1.0").split()
    info["a1"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/a2/<xmlattr>/xyz","-1.0").split()
    info["a2"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/a3/<xmlattr>/xyz","-1.0").split()
    info["a3"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))


    # Reciprocal lattice vectors
    tmp = ctx.get("Cell/b1/<xmlattr>/xyz","-1.0").split()
    info["b1"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/b2/<xmlattr>/xyz","-1.0").split()
    info["b2"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/b3/<xmlattr>/xyz","-1.0").split()
    info["b3"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))


    # K-points
    # Weights of the k-points
    tmp = ctx.get("Kmesh/weights","-1.0").split()
    weights = []
    for x in tmp:
        weights.append( float(x) )
    info["weights"] = weights

    # K-point vectors
    K = []
    tmp = ctx.get("Kmesh/k","-1.0").split()
    nk = len(tmp)/3 # the number of k-points
    for ik in xrange(nk):
        k = VECTOR(float(tmp[3*ik+0]), float(tmp[3*ik+1]), float(tmp[3*ik+2]))
        K.append(k)
    info["k"] = K       

    
    # All the bands
    all_e = []  # all bands

    for ik in range(1,info["nk"]+1):

        e_str = ctx.get("Eigenvalues/e.%i" % ik,"n").split() 
        nbnd = len(e_str)    

        norbs = len(orb_list)
        for o in orb_list:
            if o > nbnd:
                print "Orbital ", o, " is outside the range of allowed orbital indices. The maximal value is ", nbnd
                sys.exit(0)
            elif o < 1:
                print "Orbital ", o, " is outside the range of allowed orbital indices. The minimal value is ", 1
                sys.exit(0)


        e = CMATRIX(norbs,norbs)
        
        for band in orb_list:    
            mo_indx = orb_list.index(band) 
        
            ei = float(e_str[band-1])  # band starts from 1, not 0!
        
            e.set(mo_indx,mo_indx, ei * Ry2Ha, 0.0)

        all_e.append(e)

    if verbose==1:
        print " nspin = ", info["nspin"], " nk = ", info["nk"], " nbnd = ", info["nbnd"], " efermi = ", info["efermi"], \
        " alat = ", info["alat"], " omega = ", info["omega"], " tpiba = ", info["tpiba"], " tpiba2 = ", info["tpiba"]

        print " Direct lattice vectors: "
        print " a1 = ", info["a1"].x, info["a1"].y, info["a1"].z
        print " a2 = ", info["a2"].x, info["a2"].y, info["a2"].z
        print " a3 = ", info["a3"].x, info["a3"].y, info["a3"].z

        print " Reciprocal lattice vectors: "
        print " b1 = ", info["b1"].x, info["b1"].y, info["b1"].z
        print " b2 = ", info["b2"].x, info["b2"].y, info["b2"].z
        print " b3 = ", info["b3"].x, info["b3"].y, info["b3"].z

        print " K points: "
        for ik in xrange(info["nk"]):
            print ik, " weight = ", info["weights"][ik], " k = ", info["k"][ik].x, info["k"][ik].y, info["k"][ik].z

        print " Energies of the active orbitals for all k-points: "
        for ik in xrange(info["nk"]):
            print "ik = ", ik
            all_e[ik].show_matrix()
            print ""


    return info, all_e



def read_qe_wfc_info(filename, verbose=0):
    """

    This functions reads an ASCII/XML format file containing wavefunction
    and returns some descriptors

    Args:
        filename ( string ): This is the name of the file we will be reading to construct a wavefunction
        verbose ( int ): The flag controlling the amout of extra output:
        
            * 0 - no extra output (default)
            * 1 - print extra stuff
  
    Returns: 
        dictionary: a dictionary object that will contain key-value pairs with the descriptors

            * res["ngw"] ( int ): the number of plane waves needed to represent the orbital
            * res["igwx"] ( int ): the number of the G points = plane waves needed to
                represent the orbital for given k-point. Use this number when working with multiple k-points
            
            * res["nbnd"] ( int ): the number of bands (orbitals)
            * res["nspin"] ( int ): the desriptor of the spin-polarisation in the calculation 

                - 1: unpolarized
                - 2: polarized 
                - 4: non-collinear 

            * res["gamma_only"]: ( string = "T" or "F" ): T - use the Gamma-point storae trick, F otherwise
            * res["ik"] ( int ): index of the k point wfc
            * res["nk"] ( int ): the number of k-points in the wfc
        
    """   

    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()

    res = {}
    res["ngw"] = int(float(ctx.get("Info/<xmlattr>/ngw","-1.0")))     # the number of plane waves needed to represent the orbital
    res["igwx"] = int(float(ctx.get("Info/<xmlattr>/igwx","-1.0")))   # the number of the G points = plane waves needed to
                                                                      # represent the orbital for given k-point. Use this number 
                                                                      # when working with multiple k-points
    res["nbnd"] = int(float(ctx.get("Info/<xmlattr>/nbnd","-1.0")))   # the number of bands (orbitals)
    res["nspin"] = int(float(ctx.get("Info/<xmlattr>/nspin","-1.0"))) # 1 - unpolarized, 2 - polarized, 4 - non-collinear 
    res["gamma_only"] = ctx.get("Info/<xmlattr>/gamma_only","F")      # T - use the Gamma-point storae trick, T - do not use it
    res["ik"] = int(float(ctx.get("Info/<xmlattr>/ik","-1.0")))       # index of the k point wfc
    res["nk"] = int(float(ctx.get("Info/<xmlattr>/nk","-1.0")))       # the number of k-points in the wfc


    if verbose==1:
        print " ngw = ", res["ngw"], " igwx = ", res["igwx"],\
              " nbnd = ", res["nbnd"], " nspin = ", res["nspin"],\
              " gamma_only = ", res["gamma_only"],\
              " ik = ", res["ik"], " nk = ", res["nk"]

    return res




def read_qe_wfc_grid(filename, verbose=0):
    """

    This functions reads an ASCII/XML format file containing grid of G-points for given k-point
    and returns a list of VECTOR objects

    Args:
        filename ( string ): This is the name of the file we will be reading to construct a wavefunction
        verbose ( int ): The flag controlling the amout of extra output:
        
            * 0 - no extra output (default)
            * 1 - print extra stuff
  
    Returns: 
        VectorList = list of VECTOR objects: the definitions of all G points in the present calculation

    """

    ctx = Context(filename)  #("x.export/grid.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()


    # G-point vectors
    G = VECTORList() #[]

    tmp = ctx.get("grid","-1.0").split()
    ng = len(tmp)/3 # the number of G-points

    for ig in xrange(ng):
        g = VECTOR(float(tmp[3*ig+0]), float(tmp[3*ig+1]), float(tmp[3*ig+2]))
        G.append(g)

    if verbose==1:
        print "The number of G point on the grid for this k-point = ", ng
        print " G points: "
        for ig in xrange(ng):
            print ig, " g = ", G[ig].x, G[ig].y, G[ig].z
    

    return G



def read_qe_wfc(filename, orb_list, verbose=0):
    """

    This functions reads an ASCII/XML format file containing wavefunction
    and returns the coefficients of the plane waves that constitute the wavefunction

    Args:
        filename ( string ): This is the name of the file we will be reading to construct a wavefunction
        orb_list ( list of ints ): The indices of the orbitals which we want to consider. Orbitals indexing 
            at 1, not 0
        verbose ( int ): The flag controlling the amout of extra output:
        
            * 0 - no extra output (default)
            * 1 - print extra stuff
  
    Returns: 
        CMATRIX(ngw,norbs): The plane-wave expansion coefficients for all orbitals

    """


    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()

    res = read_qe_wfc_info(filename, verbose)
    ngw, nbnd, nspin, gamma_only = res["igwx"], res["nbnd"], res["nspin"], res["gamma_only"]

    if nspin==4:
        if orb_list[0] % 2 ==0:
            print "In SOC, the very first orbital index must be odd!\nExiting now..."
            sys.exit(0)
        if len(orb_list) % 2 == 1:
            print "In SOC, an even number of orbitals must be utilized!\nExiting now..."
            sys.exit(0)

    wfc_preprocess = "normalize"
    if gamma_only=="T":
        wfc_preprocess = "restore"


    norbs = len(orb_list)
    for o in orb_list:
        if o > nbnd:
            print "Orbital ", o, " is outside the range of allowed orbital indices. The maximal value is ", nbnd
            sys.exit(0)
        elif o < 1:
            print "Orbital ", o, " is outside the range of allowed orbital indices. The minimal value is ", 1
            sys.exit(0)


    coeff = CMATRIX(ngw,norbs)

    for band in orb_list:
        c = []
        all_coeff = ctx.get("Wfc."+str(band), "n").split(',')
        sz = len(all_coeff)

        for i in xrange(sz):
            a = all_coeff[i].split()
            for j in xrange(len(a)):
                c.append(a[j])
        sz = len(c)
        n = sz/2  # this should be equal to ngw

        mo_indx = orb_list.index(band) # - 1
        for i in xrange(n):
            coeff.set(i, mo_indx, float(c[2*i]), float(c[2*i+1]))



    #======== Normalize or restore (gamma-point trick) wavefunction ===============
    coeff2 = coeff #CMATRIX(ngw,norbs)

    if wfc_preprocess=="restore":

        coeff2 = CMATRIX(2*ngw-1,norbs)
        for o in xrange(norbs):
            coeff2.set(0, o, coeff.get(0, o))

            for i in xrange(1,ngw):
                coeff2.set(i, o, coeff.get(i, o))
                coeff2.set(i+ngw-1, o, coeff.get(i, o).conjugate() )

        ngw = 2*ngw - 1


    if wfc_preprocess=="normalize" or wfc_preprocess=="restore":

        if nspin==1 or nspin==2:

            for i in xrange(norbs):
                mo_i = coeff2.col(i)            
                nrm = (mo_i.H() * mo_i).get(0,0).real
                nrm = (1.0/math.sqrt(nrm))
        
                for pw in xrange(ngw):
                    coeff2.set(pw,i,nrm*coeff2.get(pw,i))

        elif nspin==4:  # spinor case

            for i in xrange(norbs/2):  # this will always be even
                mo_i_a = coeff2.col(2*i)
                mo_i_b = coeff2.col(2*i+1)

                nrm = ( (mo_i_a.H() * mo_i_a).get(0,0).real + (mo_i_b.H() * mo_i_b).get(0,0).real )
                nrm = (1.0/math.sqrt(nrm))
        
                for pw in xrange(ngw):
                    coeff2.set(pw,2*i,nrm*coeff2.get(pw,2*i))
                    coeff2.set(pw,2*i+1,nrm*coeff2.get(pw,2*i+1))
        

    return coeff2



def read_md_data(filename):
    """Read in the QE MD data stored in an XML file

    Args:
        filename (string): the name of the xml file that contains an MD data
            this function is specifically tailored for the QE output format

    Returns:
        tuple: (R, V, A, M, E), where:

        * R ( MATRIX(ndof x nsteps-1) ): coordinates of all DOFs for all mid-timesteps [Bohr]
        * V ( MATRIX(ndof x nsteps-1) ): velocities of all DOFs for all mid-timesteps [a.u. of velocity]
        * A ( MATRIX(ndof x nsteps-1) ): accelerations of all DOFs for all mid-timesteps [a.u. of acceleration]
        * M ( MATRIX(ndof x 1) ): masses of all DOFs [a.u. of mass]
        * E (list of ndof/3): atom names (elements) of all atoms 

    """

    # Default (empty) context object
    dctx = Context()

    ctx = Context(filename)      
    ctx.set_path_separator("/")
    steps = ctx.get_children("step")   
    nsteps = len(steps)
    atoms = steps[0].get_child("atomic_structure",dctx).get_child("atomic_positions", dctx).get_children("atom")
    nat = len(atoms)
    dt = float(ctx.get_child("input", dctx).get_child("ion_control", dctx).get_child("md", dctx).get_child("timestep", dctx).get("",""))
    specs = ctx.get_child("input", dctx).get_child("atomic_species", dctx).get_children("species")
    nspecs = len(specs)

    #========== Masses of elements =============
    PT = {} 
    for i in xrange(nspecs):
        name = specs[i].get("<xmlattr>/name", "X")
        mass = specs[i].get("mass", 1.0)
        PT.update({name:mass*units.amu})


    #========== Read the raw coordinates and assign masses ==========
    D = MATRIX(3*nat, nsteps) # coordinates
    f = MATRIX(3*nat, nsteps)
    M = MATRIX(3*nat, 1)
    E = []

    for t in xrange(nsteps):

        # ========== Coordinates =========
        atoms = steps[t].get_child("atomic_structure",dctx).get_child("atomic_positions", dctx).get_children("atom")

        for i in xrange(nat):        
            xyz_str = atoms[i].get("", "").split(' ')
            name = atoms[i].get("<xmlattr>/name", "X")
            D.set(3*i+0, t, float(xyz_str[0]) )
            D.set(3*i+1, t, float(xyz_str[1]) )
            D.set(3*i+2, t, float(xyz_str[2]) )

            #=========== And masses ==========
            if t==0:
                M.set(3*i+0, 0, PT[name])
                M.set(3*i+1, 0, PT[name])
                M.set(3*i+2, 0, PT[name])
                E.append(name)

        # ========== Forces  =========
        frcs = steps[t].get("forces","").split('\n')
        sz = len(frcs)
        cnt = 0
        for i in xrange(sz):        
            xyz_str = frcs[i].split()  
            if len(xyz_str)==3:
                f.set(3*cnt+0, t, float(xyz_str[0]) )
                f.set(3*cnt+1, t, float(xyz_str[1]) )
                f.set(3*cnt+2, t, float(xyz_str[2]) )
                cnt = cnt + 1                 
        
    #====== Compute velocities and coordinates at the mid-points ========
    R = MATRIX(3*nat, nsteps-1)
    V = MATRIX(3*nat, nsteps-1)
    A = MATRIX(3*nat, nsteps-1)

    for t in xrange(nsteps-1):
        for i in xrange(3*nat):    
            R.set(i, t, 0.5*(D.get(i, t+1) + D.get(i, t)) )
            V.set(i, t, (0.5/dt)*(D.get(i, t+1) - D.get(i, t)) )
            A.set(i, t, 0.5*(f.get(i, t+1) + f.get(i, t)) / M.get(i) )

    return R, V, A, M, E



def read_md_data_xyz(filename, PT, dt):
    """Read in the MD trajectory stored in an XYZ format

    Args:
        filename ( string ): the name of the xyz file that contains an MD data
        PT ( dict{ string:float } ): definition of the masses of different atomic species [masses in amu, where 1 is the mass of H]
        dt ( double ): MD nuclear integration time step [a.u. of time]

    Returns:
        tuple: (R, V, M, E), where:

        * R ( MATRIX(ndof x nsteps-1) ): coordinates of all DOFs for all mid-timesteps [Bohr]
        * V ( MATRIX(ndof x nsteps-1) ): velocities of all DOFs for all mid-timesteps [a.u. of velocity]
        * M ( MATRIX(ndof x 1) ): masses of all DOFs [a.u. of mass]
        * E (list of ndof/3): atom names (elements) of all atoms         

    """

    f = open(filename, "r")
    A = f.readlines()
    f.close()

    nlines = len(A)  # the number of lines in the file
    nat = int(float(A[0].split()[0])) # the number of atoms
    nsteps = int(nlines/(nat+2)) - 1 


    #========== Read the raw coordinates and assign masses ==========
    D = MATRIX(3*nat, nsteps) # coordinates
    M = MATRIX(3*nat, 1)
    E = []

    for t in xrange(nsteps):

        # ========== Coordinates =========
        for i in xrange(nat):        
            xyz_str = A[t*(nat+2)+2+i].split()

            name = xyz_str[0]
            D.set(3*i+0, t, float(xyz_str[1]) * units.Angst )
            D.set(3*i+1, t, float(xyz_str[2]) * units.Angst )
            D.set(3*i+2, t, float(xyz_str[3]) * units.Angst )

            #=========== And masses ==========
            if t==0:
                M.set(3*i+0, 0, PT[name]*units.amu)
                M.set(3*i+1, 0, PT[name]*units.amu)
                M.set(3*i+2, 0, PT[name]*units.amu)
                E.append(name)


    #====== Compute velocities and coordinates at the mid-points ========
    R = MATRIX(3*nat, nsteps-1)
    V = MATRIX(3*nat, nsteps-1)

    for t in xrange(nsteps-1):
        for i in xrange(3*nat):    
            R.set(i, t, 0.5*(D.get(i, t+1) + D.get(i, t)) )
            V.set(i, t, (0.5/dt)*(D.get(i, t+1) - D.get(i, t)) )

    return R, V, M, E



def read_md_data_cell(filename):
    """Read in the QE MD unit cell vectors stored in an XML file

    Args:
        filename (string): the name of the xml file that contains an MD data
            this function is specifically tailored for the QE output format

    Returns:
        MATRIX(9 x nsteps): cell coordinates for all timesteps, in the format [Bohr]: 
                 
            C(0,0) = a1.x    C(0,1) = a1.y    C(0,2) = a1.z
            C(1,0) = a2.x    C(1,1) = a2.y    C(1,2) = a2.z
            C(2,0) = a3.x    C(2,1) = a3.y    C(2,2) = a3.z

    """

    # Default (empty) context object
    dctx = Context()

    ctx = Context(filename)      
    ctx.set_path_separator("/")
    steps = ctx.get_children("step")   
    nsteps = len(steps)

    #========== Read the raw coordinates and assign masses ==========
    C = MATRIX(9, nsteps)     # cell parameters


    for t in xrange(nsteps):

        # ========== Cell =========
        cell = steps[t].get_child("atomic_structure",dctx).get_child("cell", dctx)

        a1_str = cell.get("a1", "").split(' ')
        a2_str = cell.get("a2", "").split(' ')
        a3_str = cell.get("a3", "").split(' ')

        C.set(0, t, float(a1_str[0]) );  C.set(1, t, float(a1_str[1]) );   C.set(2, t, float(a1_str[2]) );        
        C.set(3, t, float(a2_str[0]) );  C.set(4, t, float(a2_str[1]) );   C.set(5, t, float(a2_str[2]) );        
        C.set(6, t, float(a3_str[0]) );  C.set(7, t, float(a3_str[1]) );   C.set(8, t, float(a3_str[2]) );        
        
    return C




def out2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt):
    """
 
    Converts a QE output file with an MD trajectory to a bunch of input files
    for SCF calculations. These input files all have the same control settings,
    but differ in atomic coordinates
    
    Args:
        out_filename ( string ): name of the file which contains the MD trajectory
        templ_filename ( string ): name of the template file for input generation
            should not contain atomic positions! 

        prefix ( string ): the prefix of the files generated as the output
        wd ( string ): working directory where all files will be created/processed - will be created 
        t0 ( int ): defines the starting timestep to process (not all the MD timesteps may be precessed)
        tmax ( int ): defines the maximal timestep to process (not all the MD timesteps may be precessed)
        dt ( int ):  defines the spacing between frames which are written this is defined as a 
            difference between written configuration indexes. So if dt = 5, the following frames
            will be written: 0,5,10,15, etc...

    Returns:
        None: but will create a bunch of new input files in the created working directory
    """

    verbose = 0
    # Read the template file
    f_templ = open(templ_filename,"r")
    T = f_templ.readlines()
    f_templ.close()
    T_sz = len(T)

    if os.path.isdir(wd):
        pass
    else:
        # Create working directory and generate the files
        os.system("mkdir %s" % wd)

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    f_t = open("%s/tmp" % wd, "w")
    f_t.close()

    for a in f:
        #print t
#        print "is_nat=%d, nat=%d, at_line=%d, start=%d, t=%d" % (is_nat,nat,at_line,start,t)
        if start==1:
            if at_line<=nat:
                f_t.write(a)
            else:
                # A blank line just in case
                f_t.write("\n")
                f_t.close()
                t = t + 1
                start = 0
            at_line = at_line + 1

        if start==2:
            if at_line<nat:
                a_tmp = a.split()
                f_t.write(a_tmp[1]+"   "+a_tmp[6]+"   "+a_tmp[7]+"  "+a_tmp[8]+"\n")
            else:
                # A blank line just in case
                f_t.write("\n")
                f_t.close()
                t = t + 1
                start = 0
            at_line = at_line + 1


        elif start==0:
            if is_nat==0:            
                if a.find("number of atoms/cell")!=-1:
                    tmp = a.split()
                    nat = int(float(tmp[4]))
                    is_nat = 1
            else:
                if a.find("ATOMIC_POSITIONS")!=-1:
                    if t>=t0 and t<=tmax:
                        if math.fmod(t-t0,dt)==0:
                            f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                            # Write the template header
                            i = 0
                            while i<T_sz:
                                f_t.write(T[i])
                                i = i + 1
                            # Write the header for positions
                            f_t.write(a)
                            # Set flag to write positions
                            start = 1
                            at_line = 0
                        else:
                            t = t + 1
                    elif t>tmax:
                        break;
                    else:
                        t = t + 1

                elif a.find("site n.     atom                  positions (alat units)")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t>=t0 and t<=tmax:
                        if math.fmod(t-t0,dt)==0:
                            f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                            # Write the template header
                            i = 0
                            while i<T_sz:
                                f_t.write(T[i])
                                i = i + 1
                            # Write the header for positions
                            f_t.write("ATOMIC_POSITIONS (alat)\n")
                            # Set flag to write positions
                            start = 2
                            at_line = 0

                        else:
                            t = t + 1
                    elif t>tmax:
                        break;
                    else:
                        t = t + 1


    f.close()
   




def out2pdb(out_filename,T,dt,pdb_prefix):
    """

    Converts a QE output file with an MD trajectory to a bunch of pdb files
    
    Args:
        out_filename ( string ): name of the file which contains the MD trajectory.
            This file is the QE MD trajectory output
        T ( int ): defines the maximal timestep to process (not all the MD timesteps may be precessed)
        dt ( int ):  defines the spacing between frames which are written this is defined as a 
            difference between written configuration indexes. So if dt = 5, the following frames
            will be written: 0,5,10,15, etc...

    Returns:
        None: but will create a bunch of new pdb files with the trajectory info, including the unit cell params.


    Example:

        >>> QE_methods.out2pdb("x.md.out",250,25,"snaps/snap_")

       This will create MD snapshots at times 0 (input configuration), 25 (25-th nuclear configuration), 50, etc.
       the files will be collcted in the folder /snaps and named snap_0.pdb, snap_25.pdb, snap_50.pdb, etc.

    """

    dt = dt - 1
    verbose = 0

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    is_a1 = 0 # flag defining if the cell variables are set
    is_a2 = 0 # flag defining if the cell variables are set
    is_a3 = 0 # flag defining if the cell variables are set
    is_alat=0 # flag defining scaling constant for the cell
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    a1 = [0.0, 0.0, 0.0]
    a2 = [0.0, 0.0, 0.0]
    a3 = [0.0, 0.0, 0.0]

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    fr = open("tmp","w")
    fr.close()

    t = 0  # time (configuration)

    nframe = dt

    for a in f:
        if (start==1 or start==2) and t<T and nframe==dt:
            if at_line<nat:
                tmp  = a.split()
                symb = tmp[0]
                r    = [0.0, 0.0, 0.0]
                if(start==1):
                    r    = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
                    symb = tmp[0]
                elif(start==2):
                    r    = [float(tmp[6]),float(tmp[7]),float(tmp[8])]
                    symb = tmp[1]

                scl  = alat * 0.52918  # alat in a.u., coordinates will be in Angstrom
                name = ""
                elt = ""
                if len(symb)==1:
                    name = "   "+symb
                    elt = " "+symb
                elif len(symb)==2:
                    name = "  "+symb
                    elt = symb
                else:
                    name = "    "
             
                fr.write("ATOM  %5.d %s MOL M%4.d    %8.3f%8.3f%8.3f%6.2f%6.2f      XYZ1%s  \n"  % (at_line+1,name,1,scl*r[0],scl*r[1],scl*r[2],1.0,0.0,elt))
            else:
                fr.close()
                t = t + 1
                start = 0
                nframe = 0
            at_line = at_line + 1

        elif start==0:
            if is_nat==0:            
                if a.find("number of atoms/cell")!=-1:
                    tmp = a.split()
                    nat = int(float(tmp[4]))
                    #print "nat = %5d" % nat
                    is_nat = 1
            if is_a1==0:
                if a.find("a(1) =")!=-1:
                    tmp = a.split()
                    a1 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a1 = ", a1
                    is_a1 = 1
            if is_a2==0:
                if a.find("a(2) =")!=-1:
                    tmp = a.split()
                    a2 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a2 = ", a2
                    is_a2 = 1
            if is_a3==0:
                if a.find("a(3) =")!=-1:
                    tmp = a.split()
                    a3 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a3 = ", a3
                    is_a3 = 1
            if is_alat==0:
                if a.find("lattice parameter (alat)  =")!=-1:
                    tmp = a.split()
                    alat = float(tmp[4])
                    #print "alat = ", alat
                    is_alat = 1

            else:
                if a.find("ATOMIC_POSITIONS")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
                            fr = open("%s%d.pdb" % (pdb_prefix,t), "w")
                            A = math.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
                            B = math.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)
                            C = math.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)
                            AB = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]
                            AC = a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2]
                            BC = a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2]
                            gam = math.degrees(math.acos(AB/(A*B)))
                            bet = math.degrees(math.acos(AC/(A*C)))
                            alp = math.degrees(math.acos(BC/(B*C)))

                            A = A*alat*0.52918
                            B = B*alat*0.52918
                            C = C*alat*0.52918
                            fr.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (A,B,C,alp,bet,gam))

                            start = 1
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break

                elif a.find("site n.     atom                  positions (alat units)")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
                            fr = open("%s%d.pdb" % (pdb_prefix,t), "w")
                            A = math.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
                            B = math.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)
                            C = math.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)
                            AB = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]
                            AC = a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2]
                            BC = a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2]
                            gam = math.degrees(math.acos(AB/(A*B)))
                            bet = math.degrees(math.acos(AC/(A*C)))
                            alp = math.degrees(math.acos(BC/(B*C)))

                            A = A*alat*0.52918
                            B = B*alat*0.52918
                            C = C*alat*0.52918
                            fr.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (A,B,C,alp,bet,gam))

                            start = 2
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break
           
    f.close()
   



def out2xyz(out_filename,T,dt,xyz_filename):
    """
  
    This function converts the QE output file into a .xyz trajectory file. 
    No more than T steps from the out_filename file will be used

    Args: 
        out_filename ( string ): name of the file which contains the MD trajectory.
            This file is the QE MD trajectory output
        T ( int ): defines the maximal timestep to process (not all the MD timesteps may be precessed)
        dt ( int ):  defines the spacing between frames which are written this is defined as a 
            difference between written configuration indexes. So if dt = 5, the following frames
            will be written: 0,5,10,15, etc...
        xyz_filename ( string ): the name of the file to which the resulting .xyz trajectory will be written

    Returns:
        None: but will create a .xyz file with the trajectory.

    Example:

        >>> QE_methods.out2xyz("x.md.out",250,25,"snaps/traj.xyz")

        This will create the MD trajectory file in .xyz format with the snapshots takes at times 0
        (input configuration), 25 (25-th nuclear configuration), 50, etc.
        the snapshots will written in the file trahj.xyz in the folder /snaps 

    """

#    dt = dt - 1
    verbose = 0

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    is_a1 = 0 # flag defining if the cell variables are set
    is_a2 = 0 # flag defining if the cell variables are set
    is_a3 = 0 # flag defining if the cell variables are set
    is_alat=0 # flag defining scaling constant for the cell
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    a1 = [0.0, 0.0, 0.0]
    a2 = [0.0, 0.0, 0.0]
    a3 = [0.0, 0.0, 0.0]

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    fr = open(xyz_filename,"w")
    fr.close()
    fr = open(xyz_filename,"w")
    t = 0  # time (configuration)

    nframe = dt

    for a in f:
        if (start==1 or start==2) and t<T and nframe==dt:
            if at_line<nat:
                tmp  = a.split()
                symb = tmp[0]
                r    = [0.0, 0.0, 0.0]
                if(start==1):
                    r    = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
                    symb = tmp[0]
                elif(start==2):
                    r    = [float(tmp[6]),float(tmp[7]),float(tmp[8])]
                    symb = tmp[1]

                scl  = alat * 0.52918  # alat in a.u., coordinates will be in Angstrom
                fr.write("%s   %5.3f  %5.3f  %5.3f \n"  %(symb,scl*r[0],scl*r[1],scl*r[2]))
            else:
                t = t + 1
                start = 0
                nframe = 0
            at_line = at_line + 1

        elif start==0:
            if is_nat==0:            
                if a.find("number of atoms/cell")!=-1:
                    tmp = a.split()
                    nat = int(float(tmp[4]))
                    #print "nat = %5d" % nat
                    is_nat = 1
            if is_a1==0:
                if a.find("a(1) =")!=-1:
                    tmp = a.split()
                    a1 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a1 = ", a1
                    is_a1 = 1
            if is_a2==0:
                if a.find("a(2) =")!=-1:
                    tmp = a.split()
                    a2 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a2 = ", a2
                    is_a2 = 1
            if is_a3==0:
                if a.find("a(3) =")!=-1:
                    tmp = a.split()
                    a3 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a3 = ", a3
                    is_a3 = 1
            if is_alat==0:
                if a.find("lattice parameter (alat)  =")!=-1:
                    tmp = a.split()
                    alat = float(tmp[4])
                    #print "alat = ", alat
                    is_alat = 1

            else:
                if a.find("ATOMIC_POSITIONS")!=-1:
                    # Set flag to write positions
#                   print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
#                            print "Start of output"
                            fr.write("%d\n" % nat)
                            fr.write("molecule\n")

                            start = 1
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break

                elif a.find("site n.     atom                  positions (alat units)")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
#                            print "Start of output"
                            fr.write("%d\n" % nat)
                            fr.write("molecule\n")

                            start = 2
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break


    fr.close()
    f.close()
   

def xyz2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt):
    """

    Converts a xyz trajectory file with an MD trajectory to a bunch of input files
    for SCF calculations. These input files all have the same control settings,
    but differ in atomic coordinates
    
    Args:
        out_filename ( string ): name of the file which contains the MD trajectory (in xyz format)
        templ_filename ( string ): name of the template file for input generation
            should not contain atomic positions! 

        wd ( string ): working directory where all files will be created/processed - will be created 
        prefix ( string ): the prefix of the files generated as the output
        t0 ( int ): defines the starting timestep to process (not all the MD timesteps may be precessed)
        tmax ( int ): defines the maximal timestep to process (not all the MD timesteps may be precessed)
        dt ( int ):  defines the spacing between frames which are written this is defined as a 
            difference between written configuration indexes. So if dt = 5, the following frames
            will be written: 0,5,10,15, etc...

    Returns:
        None: but will create a bunch of new input files in the created working directory

    """
    verbose = 0
    # Read the template file
    f_templ = open(templ_filename,"r")
    T = f_templ.readlines()
    f_templ.close()
    T_sz = len(T)

    if os.path.isdir(wd):
        pass
    else:
        # Create working directory and generate the files
        os.system("mkdir %s" % wd)

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    f_t = open("%s/tmp" % wd, "w")
    f_t.close()

    count = 0
    # for a in the xyz file:
    for a in f:

        count += 1
        if start==0:

            b = a.strip().split()

            if is_nat==0:

                if len(b) == 1:
                    nat = int(float(b[0]))
                    is_nat = 1

            else:

                if b[0] == "Frame":

                    if t>=t0 and t<=tmax:
                        if math.fmod(t-t0,dt)==0:
                            f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                            # Write the template header
                            i = 0
                            while i<T_sz:
                                f_t.write(T[i])
                                i = i + 1
                            # Write the header for positions
                            f_t.write("ATOMIC_POSITIONS (angstrom)\n")
                            # Set flag to write positions
                            start = 1
                            at_line = 0
                        else:
                            t = t + 1
                    elif t>tmax:
                        break;
                    else:
                        t = t + 1

        elif start==1:
            if at_line<=nat-1:
                f_t.write(a)
            else:
                # A blank line just in case
                f_t.write("\n")
                f_t.close()
                t = t + 1
                start = 0
            at_line = at_line + 1

    f.close()




def get_QE_normal_modes(filename, verbosity=0):
    """

    This function reads the QE phonon calculations output files
    to get the key information for further normal modes visualization
    or other types of calculations related to normal modes  
 
    Args:  
        filename ( string ): the name of a .dyn file produced by QE code
        verbosity ( int ) to control the amount of printouts

            * 0 - no extra output (default)
            * 1 - print extra stuff
  
    Returns: 
        tuple: (Elts, R, U), where:

        * Elts ( list of nat string ): labels of all atoms, nat - is the number of atoms
        * R ( MATRIX(3*nat x 1) ): coordinates of all atoms [Angstrom]
        * U ( MATRIX(ndof x ndof) ): eigenvectors, defining the normal modes

    Example:
        >>> get_QE_normal_modes("silicon.dyn1", 1)     # verbose output
        >>> get_QE_normal_modes("Cs4SnBr6_T200.dyn1")  # not verbose output

    """

    #========= Read in the file and determine the key dimensions =======
    f   = open(filename,'r')
    A = f.readlines()
    f.close()

    tmp = A[2].split()
    nspec =  int(float(tmp[0]))  # number of types of atoms (species)
    nat = int(float(tmp[1]))     # number of atoms
    if verbosity>0:
        print "%i atoms of %i types" % (nat, nspec)


    #============= Determine the types of atoms ===============
    pfreq_indx = '(?P<freq_indx>'+rgl.INT+')'
    pAtom_type = '(?P<Atom_type>'+rgl.INT+')'+rgl.SP
    pAtom_type2 = '(?P<Atom_type2>'+rgl.INT+')'+rgl.SP
    PHRASE1 = '\'(?P<Atom_element>[a-zA-Z]+)\s+\''+rgl.SP

    last_index = 0
    E = {}
    for a in A:        
        m1 = re.search(pAtom_type + PHRASE1 + rgl.pX_val,a)
        if m1!=None:           
            ind = int(float(a[m1.start('Atom_type'):m1.end('Atom_type')]))
            elt = a[m1.start('Atom_element'):m1.end('Atom_element')]
            E.update({ind:elt})
            last_index = A.index(a)
    if verbosity>0:
        print "atom type index - element type mapping: ", E


    #============= Get the coordinates ========================
    R = MATRIX(3*nat, 1)
    Elts = []

    cnt = 0
    for i in range(last_index+1, last_index+1+nat):    
        tmp = A[i].split()
        ind = int(float(tmp[1]))
        x = float(tmp[2])
        y = float(tmp[3])
        z = float(tmp[4])

        Elts.append(E[ind])
        R.set(3*cnt+0, 0, x)
        R.set(3*cnt+1, 0, y)
        R.set(3*cnt+2, 0, z)
        cnt = cnt + 1

    if verbosity>0:
        print "Your system's elements = \n", Elts
    if verbosity>1:
        print "Your system's coordinates = \n"; R.show_matrix()


    #=========== Now look for frequencies ===============
    ndof = 3*nat
    U = MATRIX(ndof, ndof)

    pattern = 'freq \(\s+' + pfreq_indx + '\) \=\s+' + rgl.pX_val + '\[THz\] \=\s+' + rgl.pY_val + '\[cm\-1\]\s+'
    sz = len(A)
    cnt = 0
    for i in xrange(sz):        
        m1 = re.search(pattern, A[i])
        if m1!=None:
            ind = A[i][m1.start('freq_indx'):m1.end('freq_indx')] 
            freq1 = A[i][m1.start('X_val'):m1.end('X_val')] 
            freq2 = A[i][m1.start('Y_val'):m1.end('Y_val')] 
            #print ind, freq1, freq2
 
            for j in xrange(nat):
                tmp = A[i+j+1].split()
                x = float(tmp[1])
                y = float(tmp[3])
                z = float(tmp[5])

                U.set(3*j+0, cnt, x)
                U.set(3*j+1, cnt, y)
                U.set(3*j+2, cnt, z)
                
            i = i + nat
            cnt = cnt + 1

    if verbosity>1:
        print "Eigenvectors = \n"; U.show_matrix()

    return Elts, R, U
   


