#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file QE_methods.py
# This module implements functions for dealing with the outputs from QE (Quantum Espresso) package

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def cryst2cart(a1,a2,a3,r):
# auxiliary function
# convert vertor <r> in crystal (alat) coordinates with the cell defined by
# vectors a1, a2, a3, to the Cartesian coordinate xyz
    x = [0.0,0.0,0.0]
    for i in [0,1,2]:
        x[i] = a1[i]*r[0] + a2[i]*r[1] + a3[i]*r[2]

    return x



def read_qe_index(filename, orb_list, verbose=0):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the coefficients of the plane waves that constitute the
#  wavefunction
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  \param[in] upper_tag This is the name of the upper-level tag
#  \param[in] orb_list The list containing the indices of the orbitals which we want to consider
#   (indexing is starting with 1, not 0!)
#  
# Returns: a diagonal matrix (complex-valued) of orbital energies - only for those that are of interest
# in Ha = a.u. of energy


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
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the some descriptors
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  
   
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
##
#  This functions reads an ASCII/XML format file containing grid of G-points for given k-point
#  and returns a list of VECTOR objects
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#     

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
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the coefficients of the plane waves that constitute the
#  wavefunction
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  \param[in] orb_list The list containing the indices of the orbitals which we want to consider
#   (indexing is starting with 1, not 0!)
#     

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
    """
    filename (string) - the name of the xml file that contains an MD data
    this function is specifically tailored for the QE output format

    Returns:
    R ( MATRIX(ndof x nsteps-1) ) - coordinates of all DOFs for all mid-timesteps
    V ( MATRIX(ndof x nsteps-1) ) - velocities of all DOFs for all mid-timesteps
    A ( MATRIX(ndof x nsteps-1) ) - accelerations of all DOFs for all mid-timesteps
    M ( MATRIX(ndof x 1) ) - masses of all DOFs
    E (list of ndof/3) - atom names (elements) of all atoms

    All quantities are in atomic units
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
        PT.update({name:mass*amu})


    #========== Read the raw coordinates and assign masses ==========
    D = MATRIX(3*nat, nsteps)
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
            D.set(3*i+2, t, float(xyz_str[1]) )

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
                f.set(3*cnt+2, t, float(xyz_str[1]) )
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




def out2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt):
    """
    out_filename - name of the file which contains the MD trajectory
    templ_filename - name of the template file for input generation, should not contain atomic positions!
    prefix - is the prefix of the files generated at output
    wd - working directory - will be created 
    t0 and t1 - define the starting and final frames
    dt - defines the spacing between frames which are written
    this is defined as a difference between written configuration indexes:
    so if dt = 5, the following frames will be written: 0,5,10,15, etc...
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
    out_filename - name of the file which contains the MD trajectory
    this file is the QE MD trajectory output
    The function will convert this output file into pdb file (containing trajectory)
    No more than T steps from the out_filename file will be used
    dt - the difference of indexes of the frames which are written consequetively
    such that if you dt = 5 it will write frames 0,5,10,15,etc. with 0 - being the input configuration

    Example of usage:
    > out2pdb.convert("x.md.out",250,25,"snaps/snap_")
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
    out_filename - name of the file which contains the MD trajectory
    this file is the QE MD trajectory output
    The function will convert this output file into xyz file (containing trajectory)
    No more than T steps from the out_filename file will be used
    dt - number of steps between output frames, so dt = 5 will output frames 0, 5, 10, 15, etc.
    xyz_filename - is the prefix of the file to which the result is written

    Example of usage:
    > out2xyz.convert("x.md.out",250,25,"snaps/traj.xyz")
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
    out_filename - name of the file which contains the xyz (MD) trajectory
    templ_filename - name of the template file for input generation, should not contain atomic positions!
    prefix - is the prefix of the files generated at output
    wd - working directory - will be created 
    t0 and t1 - define the starting and final frames
    dt - defines the spacing between frames which are written
    this is defined as a difference between written configuration indexes:
    so if dt = 5, the following frames will be written: 0,5,10,15, etc...
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

