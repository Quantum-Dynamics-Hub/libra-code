#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: fix_motion
   :platform: Unix, Windows
   :synopsis: This module implements functions for removing translation of CM and
       rotation around it

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

from . import units
from . import build
import libra_py.packages.qe.methods as QE_methods
from . import normal_modes
from . import LoadPT
from . import LoadMolecule
from . import nve_md


def permutation1(case):
    """
    This function returns one of the 6 permutation matrices that
    describe possible orderings/directions of the 3 coordinate axes.

    Args:
        case ( int ): selector of the permutation

    Returns:
        MATRIX3x3: the matrix that encodes the permutation

    """
    
    U = MATRIX3x3()
        
    if case==0:
        # U[0] * diag(I1,I2,I3) * U[0] = diag(I1,I2,I3)
        # det = 1
        # identity
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;
        
    elif case==1:
        # U[1] * diag(I1,I2,I3) * U[1] = diag(I2,I1,I3)
        # det = 1 
        # x -> y', y -> x', z -> -z'
        U.xx = 0.0; U.xy = 1.0; U.xz = 0.0;
        U.yx = 1.0; U.yy = 0.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;

    elif case==2:
        # U[2] * diag(I1,I2,I3) * U[2] = diag(I1,I3,I2)
        # det = 1 
        # x -> -x', y -> z', z -> y'
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 0.0; U.yz = 1.0;
        U.zx = 0.0; U.zy = 1.0; U.zz = 0.0;

    elif case==3:
        # U[3] * diag(I1,I2,I3) * U[4] = diag(I3,I1,I2)
        # det = 1
        # x -> z', y -> x', z -> y'
        U.xx = 0.0; U.xy = 0.0; U.xz = 1.0;
        U.yx = 1.0; U.yy = 0.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 1.0; U.zz = 0.0;

    elif case==4:
        # U[4] * diag(I1,I2,I3) * U[3] = diag(I2,I3,I1)
        # det = 1
        # x -> y', y -> z', z -> x'
        U.xx = 0.0; U.xy = 1.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 0.0; U.yz = 1.0;
        U.zx = 1.0; U.zy = 0.0; U.zz = 0.0;

    elif case==5: 
        # U[5] * diag(I1,I2,I3) * U[5] = diag(I3,I2,I1)
        # det = 1 
        # x -> z', y -> -y', z -> x'
        U.xx = 0.0; U.xy = 0.0; U.xz = 1.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 1.0; U.zy = 0.0; U.zz = 0.0;
    
    return U


def permutation(case):
    """
    This function returns one of the 6 permutation matrices that
    describe possible orderings/directions of the 3 coordinate axes.

    Args:
        case ( int ): selector of the permutation

    Returns:
        MATRIX3x3: the matrix that encodes the permutation

    """
    
    U = MATRIX3x3()
        
    if case==0:
        # x -> x    y -> y      z -> z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;
        
    elif case==1:
        # x <-> y     z -> z
        U.xx = 0.0; U.xy = 1.0; U.xz = 0.0;
        U.yx = 1.0; U.yy = 0.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;

    elif case==2:
        # x <-> z     y -> y
        U.xx = 0.0; U.xy = 0.0; U.xz = 1.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 1.0; U.zy = 0.0; U.zz = 0.0;

    elif case==3:
        # y <-> z     x -> x
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 0.0; U.yz = 1.0;
        U.zx = 0.0; U.zy = 1.0; U.zz = 0.0;

    
    return U


def permutation2(case):
    """
    This function returns one of the 6 permutation matrices that
    describe possible orderings/directions of the 3 coordinate axes.

    Args:
        case ( int ): selector of the permutation

    Returns:
        MATRIX3x3: the matrix that encodes the permutation

    """
    
    U = MATRIX3x3()
        
    if case==0:
        # x -> x    y -> y      z -> z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;
        
    elif case==1:
        # x <-> y     z -> -z
        U.xx = 0.0; U.xy = 1.0; U.xz = 0.0;
        U.yx = 1.0; U.yy = 0.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;

    elif case==2:
        # x <-> z     y -> -y
        U.xx = 0.0; U.xy = 0.0; U.xz = 1.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 1.0; U.zy = 0.0; U.zz = 0.0;

    elif case==3:
        # y <-> z     x -> -x
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 0.0; U.yz = 1.0;
        U.zx = 0.0; U.zy = 1.0; U.zz = 0.0;

    elif case==4:
        # x->y, y->z, z->x
        U.xx = 0.0; U.xy = 0.0; U.xz = 1.0;
        U.yx = 1.0; U.yy = 0.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 1.0; U.zz = 0.0;

    elif case==5:
        # x->z, z->y, y->x
        U.xx = 0.0; U.xy = 1.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 0.0; U.yz = 1.0;
        U.zx = 1.0; U.zy = 0.0; U.zz = 0.0;


    
    return U



def reflection(case):

    U = MATRIX3x3()
        
    if case==0:
        # x - > x, y -> y, z -> z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;


    elif case==1:
        # x - > -x, y -> y, z -> z
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;

    elif case==2:
        # x - > x, y -> -y, z -> z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;

    elif case==3:
        # x - > x, y -> y, z -> -z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;


    elif case==4:
        # x - > -x, y -> -y, z -> z
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz = 1.0;

    elif case==5:
        # x - > -x, y -> y, z -> -z
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy = 1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;

    elif case==6:
        # x - > x, y -> -y, z -> -z
        U.xx = 1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;


    elif case==7:
        # x - > -x, y -> -y, z -> -z
        U.xx =-1.0; U.xy = 0.0; U.xz = 0.0;
        U.yx = 0.0; U.yy =-1.0; U.yz = 0.0;
        U.zx = 0.0; U.zy = 0.0; U.zz =-1.0;

    return U

       

def measure(S0, S1):
    """
    Compute the measure of the two systems' dis-orientation
    It is zero if the two systems are ideally parallel to each other

    Args:
        S0 ( System ): first system
        S1 ( System ): second system
        
    Returns:
        ( double ): measure of two systems' dis-orientation

    """
    nat = S0.Number_of_atoms
    res, dot = 0.0, 0.0
    for n in range(0,nat):
        
        dr0 = S0.Atoms[n].Atom_RB.rb_cm -  S0.Fragments[0].Group_RB.rb_cm 
        dr1 = S1.Atoms[n].Atom_RB.rb_cm -  S1.Fragments[0].Group_RB.rb_cm 
        
        dr0.normalize()
        dr1.normalize()
        l = cross(1.0, dr0, dr1)

        res = res + l.length2()        
        dot = dot + dr0 * dr1

    res = res / float(nat)
    dot = dot / float(nat)
    
    return res,  dot


def measure2(S0, S1):
    """
    Compute the measure of the two systems' dis-orientation
    It is zero if the two systems are ideally parallel to each other

    Args:
        S0 ( System ): first system
        S1 ( System ): second system
        
    Returns:
        ( double ): measure of two systems' dis-orientation

    """
    nat = S0.Number_of_atoms
    res, dot = 0.0, 0.0
    for n in range(0,nat):
        
        dr = S0.Atoms[n].Atom_RB.rb_cm - S1.Atoms[n].Atom_RB.rb_cm       
        dot = dot + dr.length2()


        dr0 = S0.Atoms[n].Atom_RB.rb_cm -  S0.Fragments[0].Group_RB.rb_cm 
        dr1 = S1.Atoms[n].Atom_RB.rb_cm -  S1.Fragments[0].Group_RB.rb_cm         
        dr0.normalize()
        dr1.normalize()
        res = res + dr0 * dr1

    res = res / float(nat)
    dot = dot / float(nat)
     
    return dot, res



def remove_translation(S0, S1):
    """
    Translates the system ```S1``` such that it's center of mass 
    coincides with that of the system ```S0```

    Args:
        S0 ( System ): first system
        S1 ( System ): second system
        
    Returns:
        None: but the system ```S1``` is changed

    """


    R0 = S0.Fragments[0].Group_RB.rb_cm 
    R1 = S1.Fragments[0].Group_RB.rb_cm 
    dR = R0 - R1
    dr = dR.length()

    #print dR.x, dR.y, dR.z
    S1.TRANSLATE_FRAGMENT(dr, dR, 1)


def remove_rotation(S0, S1, verbose=0):
    """
    Rotates the system ```S1``` such that it is oriented as similar to 
    the system ```S0``` as possible

    Args:
        S0 ( System ): first system
        S1 ( System ): second system
        verbose ( int ): the level of debug info printout [ default: 0 - no printout ]
        
    Returns:
        None: but the system ```S1``` is changed

    """
    I = MATRIX3x3(); I.identity()

    Q0 = QUATERNION()
    Q1 = QUATERNION()
    Q2 = QUATERNION()
    U = MATRIX3x3()

    MATRIX_TO_QUATERNION(S0.Fragments[0].Group_RB.rb_A_I_to_e, Q0);
    MATRIX_TO_QUATERNION(S1.Fragments[0].Group_RB.rb_A_I_to_e, Q1);

    if verbose > 0:
        print("Q0 = ", Q0.Lt, Q0.Lx, Q0.Ly, Q0.Lz)
        print("Q1 = ", Q1.Lt, Q1.Lx, Q1.Ly, Q1.Lz)

    # If one of the bodies is totally-symmetric, don't do anything
    # bacause we can't e sure about the directions
    if Q0.vect().length2() * Q1.vect().length2() > 1e-15:
         
        if verbose > 0:       
            print("Measure before rotation = ", measure2(S0, S1) )

        # Lets consider all possible axes permutations and compute the
        # alignment measure for each one
        mes = []
        indices = []

        for i in range(0,8):            # reflections
            for j in range(0,6):        # permutations
#        for i in [0]:            # reflections
#            for j in [0]:        # permutations


                #Q2 = Q0 * Q1.inverse()
                #QUATERNION_TO_MATRIX(Q2, U);


                P = reflection(i) * permutation2(j)  # reflection(i) *
                #P = permutation2(j) * reflection(i)
                A1 = S1.Fragments[0].Group_RB.rb_A_I_to_e * P    
                U = S0.Fragments[0].Group_RB.rb_A_I_to_e * A1.T() #A1.inverse()

                """
                dU = A1
                print "%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz)
                print "%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz)
                print "%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz)

                dU = S0.Fragments[0].Group_RB.rb_A_I_to_e
                print "%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz)
                print "%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz)
                print "%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz)

                dU = S0.Fragments[0].Group_RB.rb_A_I_to_e * S0.Fragments[0].Group_RB.rb_A_I_to_e_T
                print "%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz)
                print "%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz)
                print "%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz)

                dU = A1 * A1.T()
                print "%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz)
                print "%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz)
                print "%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz)

                dU = (U*U.T() - I)
                print "%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz)
                print "%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz)
                print "%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz)

                print "%8.5f  %8.5f  %8.5f " % (U.xx, U.xy, U.xz)
                print "%8.5f  %8.5f  %8.5f " % (U.yx, U.yy, U.yz)
                print "%8.5f  %8.5f  %8.5f " % (U.zx, U.zy, U.zz)
                """


                if math.fabs(U.Determinant()) < 1e-10:
                    print("WARNING: 1) Problems with inversion: ", U.Determinant() )
                    sys.exit(0)

                dU = (U*U.T() - I)
                if (dU * dU.T()).tr() > 1e-10:
                    print("WARNING: 2) Problems with inversion: ", (dU * dU.T()).tr())
                    print("%8.5f  %8.5f  %8.5f " % (dU.xx, dU.xy, dU.xz))
                    print("%8.5f  %8.5f  %8.5f " % (dU.yx, dU.yy, dU.yz))
                    print("%8.5f  %8.5f  %8.5f " % (dU.zx, dU.zy, dU.zz))

                    sys.exit(0)

                s = System(S1)
                s.ROTATE_FRAGMENT(U, 1)
    
                if verbose > 1:
                    print("%8.5f  %8.5f  %8.5f " % (U.xx, U.xy, U.xz) )
                    print("%8.5f  %8.5f  %8.5f " % (U.yx, U.yy, U.yz) )
                    print("%8.5f  %8.5f  %8.5f " % (U.zx, U.zy, U.zz) )
                
                mes.append( measure2(S0, s) )
                indices.append( [i, j] )

                if verbose > 0:
                    print( "Measure after reflection %i and permutaiton %i = " % (i, j), mes[6*i+j] )


                # Rotate the system S1 back to its original configuration
                # before applying the next ordering 
                #S1.ROTATE_FRAGMENT(U.inverse(), 1)

        # Determine the permutation for which the alignment of the
        # two configurations is the best
        mes_indx, min_dist, max_dot = 0, mes[0][0], mes[0][1]
        sz = len(mes)
        for i in range(1,sz):
            if mes[i][0] < min_dist and mes[i][1] > max_dot:
                min_dist = mes[i][0]
                max_dot = mes[i][1]
                mes_indx = i


        if verbose > 0:
            print("Selected transformation = ", mes_indx, indices[mes_indx] )

        #QUATERNION_TO_MATRIX(Q2, U);

        # Apply the rotation with the correct permutation    
        P = reflection(indices[mes_indx][0]) * permutation2(indices[mes_indx][1])   
        #P = permutation2(indices[mes_indx][1]) * reflection(indices[mes_indx][0])    
        A1 = S1.Fragments[0].Group_RB.rb_A_I_to_e * P    
        U = S0.Fragments[0].Group_RB.rb_A_I_to_e * A1.T()  #A1.inverse()


        S1.ROTATE_FRAGMENT(U, 1)

        if verbose > 0:
            print("Measure after rotation = ", measure2(S0, S1) )
            MATRIX_TO_QUATERNION(S1.Fragments[0].Group_RB.rb_A_I_to_e, Q1);
            print("Q1(updated) = ", Q1.Lt, Q1.Lx, Q1.Ly, Q1.Lz)


def process_xyz(filename, PT, itime, ftime, verbose=0):
    """
    This function will read xyz file to create the coordinates for a range of 
    configurations and it will remove the translation of COM and rotation of 
    the configuration. 

    Args:
        filename ( string ): the name of the xyz file to read
        PT ( dictionary ): periodic table, which sets up the atomic masses [ in Daltons ]
            for all elements that one meets in the xyz file
        itime ( int ): the initial timeframe to process
        ftime ( int ): the final timefram to process

    Returns:
        None

    """

    # Read in the xyz file 
    R, E = QE_methods.read_md_data_xyz2(filename, PT) # R  in a.u. 

    # Determine dimensions
    ndof = R.num_of_rows
    nat = int(ndof/3)
    nsteps = R.num_of_cols
    if verbose>0:
        print("nsteps = ", nsteps)
        print("ndos = ", ndof)

    # Setup masses
    M = MATRIX(ndof, 1)
    for at in range(0,nat):    
        M.set(3*at+0, 0, PT[E[at]])
        M.set(3*at+1, 0, PT[E[at]])
        M.set(3*at+2, 0, PT[E[at]])

    # Setup the universe    
    U = Universe()
    LoadPT.Load_PT(U, "elements.dat", 0)

    # Build chemical systems corresponding to each geometry 
    # also find the reference system - the one for which the directions are 
    # not of zero length
    indx, found = 0, False
    Systems = []
    for frame in range(itime, ftime+1):
        # Build systems
        s = build.make_system(R, E, frame, U)        
        Systems.append(s)

        # Check if it is a good reference
        if found == False:
            Q = QUATERNION()
            MATRIX_TO_QUATERNION(s.Fragments[0].Group_RB.rb_A_I_to_e, Q)
            if Q.vect().length2() > 1e-15:
                indx = frame - itime
                found = True

        # Optionally, print out some info           
        if verbose>2:
            s.Fragments[0].Group_RB.show_info()

    # Now, we are ready to process all the systems and convert them to xyz
    xyz0 = ""
    xyz = ""
    nframes = len(Systems)
    for frame in range(0,nframes):

        # To xyz - before
        xyz1 = nve_md.syst2xyz(Systems[frame])
        xyz0 = xyz0 + xyz1

        if frame!=indx:
            remove_translation(Systems[indx], Systems[frame])        
            remove_rotation(Systems[indx], Systems[frame], verbose)

        # To xyz - after
        xyz1 = nve_md.syst2xyz(Systems[frame])
        xyz = xyz + xyz1


    # Convert the Systems into R
    R  = MATRIX(ndof, nframes)
    for frame in range(0,nframes):    

        # To MATRIX
        for at in range(0,nat):    
            R.set(3*at+0, frame, Systems[frame].Atoms[at].Atom_RB.rb_cm.x)
            R.set(3*at+1, frame, Systems[frame].Atoms[at].Atom_RB.rb_cm.y)
            R.set(3*at+2, frame, Systems[frame].Atoms[at].Atom_RB.rb_cm.z)
       
    return R, xyz0, xyz

