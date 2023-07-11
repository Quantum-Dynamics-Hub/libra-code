#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  This code runs a lattice Monte Carlo
#
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


def make_doped(nat, fract, base_elt, dope_elt):
    """
    nat - integer
    fract - float
    base_elt - integer
    dope_elt - integer

    Creates a system of <nat> sites occupied by the <base_elt>
    a <fract> fraction ([0,1]) of these sites is then replaced
    by the doping element <dope_elt>

    Return: a list of labels <L> showing the occupations of all the sites
    and a list of indices, <dope_sites>, of the sites occupied by the dopants 
    """

    nnew = int(fract*nat)

    dope_sites = intList()
    randperm(nnew, nat, dope_sites)

    L = []
    for i in xrange(nat):
        L.append(base_elt)

    for i in xrange(nnew):
        L[ dope_sites[i] ] = dope_elt
         
 
    return L



def energy(connect, sites, H1, H2):
    """ 
    connect - list of [int, [list of ints] ] elements
    sites - list of integers that specify which species type occupy each site
    H1 - MATRIX containing parameters of on-site energy contribution
    H2 - MATRIX containing parameters of interactions of different types

    This function computes the energy of a lattice occupied by the sites <L>        
    and connected according to <connect> object
          
    """
   
    en = 0.0
    sz = len(connect)

    #===== 1-body contributions ======
    for i in xrange(sz):
        en = en + H1.get(sites[i],sites[i])
    
    #===== 2-body contributions ======

    for i in xrange(sz):
        neib = connect[i][1]

        for j in neib:
            en = en + 0.5 * H2.get(sites[i],sites[j])
            
    return en



def swap_sites(n, sites, dopant_indx):
    """
    This function will attempt to swap <n> sites occupied by
    a dopant with <dopant_indx> index with the <n> sites occupied
    by any other type of specied
    <sites> - is a list of integers, each integer corresponds to the
    type of atom that occupies each site
    """

    nat = len(sites)

    # Determine the doped and undoped sites:
    doped_sites = []
    undoped_sites = []
    for i in xrange(nat):
        if sites[i]==dopant_indx:
           doped_sites.append(i)
        else:
           undoped_sites.append(i)


    ndop = len(doped_sites)

    # Now many pairs to swap
    npic = n
    if npic>ndop:
        npic = ndop
        if npic>n-ndop:
            npic = n-ndop

    # Pick doped and undoped sites
    pic_doped = intList()
    randperm(npic, ndop, pic_doped)

    pic_undoped = intList()
    randperm(npic, nat-ndop, pic_undoped)

    # Make a copy
    new_sites = list(sites)    

    # Swap these npic pairs of sites:
    for n in xrange(npic):
        i, j = doped_sites[ pic_doped[n] ],  undoped_sites[ pic_undoped[n] ]
        si = new_sites[i];  new_sites[i] = new_sites[j];  new_sites[j] = si;

    return new_sites




def mc_move(connect, sites, dopant_indx, n_move_at, H1, H2, oldE, kT, rnd):
    """
    This function attempts to perform a MC move of <n_move_at> atoms - we only consider
    swapping the sites occupied by a dopant with index <dopant_indx> in the global array
    of sites <sites>. We also need a connectivity object <connect> which lists neighbors of 
    each atom. <H1> and <H2> are the matrices that specify the on-site and site-site interaction
    energies, respectively. <oldE> is the energy of the previous configuration, from which we
    attempt to hop. <kT> is an effective temperature, and <rnd> is an object of needed for
    random number generation.
    
    """

    new_sites = swap_sites(n_move_at, sites, dopant_indx)

    newE = energy(connect, new_sites, H1, H2)

    acc_rat = 1.0 
    argg = (newE - oldE)/kT
    if argg>0:
        if argg<50:
            acc_rat = math.exp(-(newE - oldE)/kT)
        else:
            acc_rat = 0.0

    ksi = rnd.uniform(0.0, 1.0)
  
    if ksi<acc_rat:
        # accept the move
        return new_sites, newE

    else:
        return sites, oldE 




def main():

    #==========================================================================    
    #             Prepare system and initialize essential variables
    #==========================================================================

    Angst = 1.889725989     
    prefix = "ar"

    # xyz file contains coordinates in Angstrom, but 
    # they are converted to Bohrs when read
    L0, R0 = build.read_xyz(prefix+".xyz")  


    Nx, Ny, Nz = 30, 30, 1
    a = VECTOR(1.0,   0.0,   0.0) * Angst # convert to a.u.
    b = VECTOR(0.0,   1.0,   0.0) * Angst # convert to a.u.
    c = VECTOR(0.0,   0.0,   1.0) * Angst # convert to a.u.
    L1, R1 = build.generate_replicas_xyz2(L0, R0, a, b, c, Nx, Ny, Nz)
    
    
    masses = []
    MaxCoord = []
    for i in xrange(len(L1)):
        masses.append(1.0)
        MaxCoord.append(6)
    
    # Compute the connectivities   
    Rcut = 1.1 * Angst
    tv1 = a * Nx
    tv2 = b * Ny
    tv3 = c * Ny
    pbc_opt = "ab"
    connect, line, unsorted_pairs_out = autoconnect.autoconnect_pbc(R1, MaxCoord, Rcut, tv1, tv2, tv3, pbc_opt)
    
    
    # Create an empty chemical system object
    syst = System()
    build.add_atoms_to_system(syst, L1, R1, VECTOR(1.0, 1.0, 1.0), masses, "elements.dat")
    


    #==========================================================================    
    #                               Run MC
    #==========================================================================
    # Simulation parameters
    nat = syst.Number_of_atoms
    nsnaps = 100
    nsteps = 500
    kT = 0.01
    rnd = Random()
    fract = 0.2
    dopant_indx = 1
    n_move_at = 1
    mapping = ["Ar", "Br"]

    H1 = MATRIX(2,2)
    H1.set(0,0, 0.0)
    H1.set(1,1, 0.0)
    
    H2 = MATRIX(2,2)
    H2.set(0,0,  0.00);  H2.set(0,1, -0.05)
    H2.set(1,0, -0.05);  H2.set(1,1, -0.20)

 
    # Create the dopped lattice
    sites = make_doped(nat, fract, 0, dopant_indx)


    en = energy(connect, sites, H1, H2)
    print "initial energy = ", en

    snap = 0
    for i in xrange(nat):
        syst.Atoms[i].Atom_element = mapping[ sites[i] ] 
    syst.print_ent("2D/step%i" % snap)


    # Continue the MC moves and print out the lattices after every nsteps
    f = open("2D/energy.txt", "w")
    f.write("%5i  %8.5f \n" % (snap, en))
    f.close()


    for snap in xrange(1,nsnaps):
        for step in xrange(nsteps):

            sites, en = mc_move(connect, sites, dopant_indx, n_move_at, H1, H2, en, kT, rnd)

        for i in xrange(nat):
            syst.Atoms[i].Atom_element = mapping[ sites[i] ] 
        
        syst.print_ent("2D/step%i" % snap)

        print "step= ", snap, " energy= ", en
        f = open("2D/energy.txt", "a")
        f.write("%5i  %8.5f \n" % (snap, en))
        f.close()



main()    

