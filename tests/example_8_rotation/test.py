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

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def shift(rb, params):
    """ 
    Returns a shift vector needed to translate the rb's cm to 
    a particular location, so to define the rotation center
    """


    if params["pivot_option"]=="atom":
        """ 
        Rotate around specified atom (center)
        """

        indx = params["pivot_atom"]
        shft = rb.get_center_in_global_frame(indx) - rb.rb_cm

    elif params["pivot_option"]=="arbitrary":
        """ 
        Rotate around an arbitrary pivot point
        """

        ppt = params["pivot_point"]
        shft = ppt - rb.rb_cm

    return shft



def main(centers, masses, params, filename):
    """
    This function illustrates a rotation of the system of material points 
    evolving as a single rigid body. 
    Parameters:
    centers [list of VECTOR objects] represent the coordinates of all points
    masses [list of double] represent the masses of all points
    filename ["anything-here.xyz"] the file where the produced trajectory will be stored
    params [dictionary] - defines simulation protocol (parameters)
    """

     
    """ 
    Create and initialize a rigid body with a given number of points and distribution of masses
    """
    rb = RigidBody()
    N = len(centers)
    rb.init(N,masses,centers)
    rb.show_info()


    lab_dir = VECTOR(0.0, 0.0, 1.0)   
    body_dir = VECTOR(0.0, 0.0, 1.0)   


    """ 
    Use the system object just for the visualization purposes
    """
    U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")
    syst = System()
    for n in xrange(N):
        syst.CREATE_ATOM( Atom(U,  {"Atom_element":"C","Atom_cm_x":0.0,"Atom_cm_y":0.0,"Atom_cm_z":0.0})  )
        syst.Atoms[n].Atom_RB.rb_cm = rb.get_center_in_global_frame(n)


    """ 
    Erase the content of the output file, if it wasn't empty 
    and print the first (initial) geometry
    """
    f = open(filename, "w")
    f.close()    
    syst.print_xyz(filename,0)


    """ 
    Do the rotation and print out geometries
    """
    for i in xrange(1, params["nsteps"]):

        rb.shift_position(shift(rb, params))

        if params["rot_type"]=="lab":
            rb.Rotate_I(params["rot_amt"], params["rot_dir"]) 
        elif params["rot_type"]=="body":
            rb.Rotate_e(params["rot_amt"], params["rot_dir"]) 

        rb.shift_position(-1.0*shift(rb, params))

        for n in xrange(N):
            syst.Atoms[n].Atom_RB.rb_cm = rb.get_center_in_global_frame(n)    
        syst.print_xyz(filename,i)


#=========== Run all the tests ===============
# Rotations around atoms:
rot_dirs = [VECTOR(1.0, 0.0, 0.0), VECTOR(0.0, 1.0, 0.0), VECTOR(0.0, 0.0, 1.0), VECTOR(1.0, 1.0, 0.0)]

count = 0
for atom in [1,2,3,4]:
    for rd in [0,1,2]: #,3]:
        for rot_type in ["lab"]: #, "body"]:

            main([VECTOR(0.0,0.0,0.0),  VECTOR(2.0,0.0,0.0), VECTOR(0.0,2.0,0.0), VECTOR(0.0,0.0,2.0)], 
                 [1.0, 1.0, 1.0, 1.0],
                 {"nsteps":100, "pivot_option":"atom", "pivot_atom":atom, "rot_amt":0.1, "rot_type":rot_type, "rot_dir":rot_dirs[rd]},
                 "case-%i-atom%i-dir%i-%s.xyz" % (count, atom, rd, rot_type )
                )

            count = count + 1

