import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import LoadPT  # for Load_PT

def init_system(label, R, g, rnd, T, sigma, df, data_file="elements.dat"):
##
# This function creates and instance of the System object which will be used in classical MD.
#
# \param[in] label The list of atom names (list of strings)
# \param[in] R     The list of coordinates of all atoms (list of VECTORs), units: Bohr
# \param[in] g     The list of gradients on all atoms (list of VECTORs), units: Ha/Bohr
# \param[in] rnd   Random number generator object
# \param[in] T     Target temperature used to initialize momenta of atoms. units: K
# \param[in] sigma The magnitude of a random displacement of each atom from its center, units: Bohr
# \param[in] df    Controls additional (debug) printing
# \param[in] data_file The file containing information about elements.
#

    # Create Universe and populate it
    U = Universe();   LoadPT.Load_PT(U, data_file, 0)
    syst = System()

    sz = len(label)
    for i in xrange(sz):
        atom_dict = {}
        atom_dict["Atom_element"] = label[i]
        atom_dict["Atom_cm_x"] = R[i].x + sigma*rnd.normal()
        atom_dict["Atom_cm_y"] = R[i].y + sigma*rnd.normal()
        atom_dict["Atom_cm_z"] = R[i].z + sigma*rnd.normal()

        if df:
            print "CREATE_ATOM ",atom_dict["Atom_element"]

        at = Atom(U, atom_dict)
        at.Atom_RB.rb_force = VECTOR(-g[i].x, -g[i].y, -g[i].z)


        syst.CREATE_ATOM(at)

    if df:
        print "Number of atoms in the system = ", syst.Number_of_atoms
        syst.show_atoms()


    # Initialize random velocity at T(K) using normal distribution
    syst.init_atom_velocities(T)

    return syst

