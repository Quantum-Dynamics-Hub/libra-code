#*********************************************************************************
#* Copyright (C) 2019-2020 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: scan
   :platform: Unix, Windows
   :synopsis: This module implements functions for 
       performing various types of calculations related to PES scans  

.. moduleauthor:: Alexey V. Akimov 
  
"""


import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import units


def coords2xyz(labels, q, itraj, unit_conversion_factor=1.0):
    """

    This function creates a string containing an xyz-formatted atomistic information.

    Args:
        labels ( list of strings ): the labels of the atoms of the system, the labels are in a 
            specific order that corresponds to the atomic info in an xyz file. The length of this 
            list is `nat` - the number of atoms.
        q ( MATRIX(ndof, ntraj) ): the coordinates of all DOFs for all trajectories, here `ndof` = 3 `natoms` 
        itraj ( int ): the index of the trajectory that we are interested in.
        unit_conversion_factor ( double ): the conversion factor to convert input in `q` unit to the new units
            expected in the xyz file. In a common situation, when `q` uses Bohr units, and xyz is expected to 
            be in Angstroms, use the `unit_conversion_factor = 1.0/units.Angst` 
 
    Returns:
        string: the string representation of the xyz file that is made of the provided input geomtry/atomic labels

    """

    natoms = int(q.num_of_rows / 3)
        
    res = ""
    for iatom in range(natoms):
        label = labels[iatom]
        x = unit_conversion_factor * q.get(3*iatom+0, itraj)
        y = unit_conversion_factor * q.get(3*iatom+1, itraj)
        z = unit_conversion_factor * q.get(3*iatom+2, itraj)        
        res = res + F"{label}  {x} {y} {z}\n"
        
    return res


def make_path_xyz2(labels, path, itraj, unit_conversion_factor=1.0):
    """
    This function creates a string containing an xyz-formatted trajectory (multiple geometries)

    Args:
        labels ( list of strings ): the labels of the atoms of the system, the labels are in a 
            specific order that corresponds to the atomic info in an xyz file. The length of this 
            list is `nat` - the number of atoms.
        path ( list of nsteps MATRIX(ndof, ntraj) objects ): the coordinates of all DOFs for all trajectories for `nsteps` frames
            (geometries), here `ndof` = 3 `natoms` 
        itraj ( int ): the index of the trajectory that we are interested in.
        unit_conversion_factor ( double ): the conversion factor to convert input in `q` unit to the new units
            expected in the xyz file. In a common situation, when `q` uses Bohr units, and xyz is expected to 
            be in Angstroms, use the `unit_conversion_factor = 1.0/units.Angst` 
 
    Returns:
        string: the string representation of the xyz file that is made of the provided input geomtry/atomic labels 
    
    """

    nsnaps = len(path)
    natoms = len(labels)

    md_xyz = ""
    for isnap in range(nsnaps):
        snap = coords2xyz(labels, path[isnap], itraj, unit_conversion_factor)    
        md_xyz = md_xyz + F"{natoms}\nsnapshot {isnap}\n{snap}"
    
    return md_xyz



def make_path_xyz(R0, R1, E, s0=0.0, s1=1.0, npts=2, S0=0.0, S1=1.0):
    """
 
    This function generates an xyz file with a "path"
    connecting one geometry to another (like in NEB), but
    also extends the range of the scan beyond the initial
    limiting points, allowing the put extra points along the
    scan direction

    Args:
        R0 ( MATRIX(3*nat, 1) ): coordinates of all DOFs for initial input geometry [ Bohr ]
        R1 ( MATRIX(3*nat, 1) ): coordinates of all DOFs for final input geometry [ Bohr ]
        E (list of nat strings): atom names (elements) of all atoms         
        s0 ( double ): mapping of the R0 point onto the line of search [ default: 0.0 ]
        s1 ( double ): mapping of the R1 point onto the line of search [ default: 1.0 ]
        npts ( int ): the number of images (including the boundary points) we would have between the
            input points. This number defines the spacing between the images on the path.
        S0 ( double ): mapping of the actual initial point onto the line of search [ default: 0.0 ]
        S1 ( double ): mapping of the actual final point onto the line of search [ default: 1.0 ]

    Returns:
        tuple: ( R, xyz, s_axis), where:

            * R ( MATRIX(3*nat, Npts) ): coordinates of all points on the scan path [ Bohr ]
            * xyz ( string ): the xyz file content with all the points [ Angstrom ]
            * s_axis ( list of Npts doubles): values of the "reaction" coordinate
        
    """

    # Sanity check:
    if s1<=s0:
        print("Error: s0 (%5.3f) should be larger than s1(%5.3f)! Exiting...\n" % (s0, s1))
        sys.exit(0)

    # Determine the number of points in the the actual path:
    ndof = R0.num_of_rows
    nat = int(ndof / 3)

    print("ndof = ", ndof)
    print("nat = ", nat)

    ds = (s1-s0)/float(npts)
    Npts = int(abs((S1 - S0))/ds) + 1
    print("ds = ", ds)
    print("Npts = ", Npts)
    
    # Allocate memory
    R = MATRIX(ndof, Npts)

    s_axis = []

    # Compute the mapping
    for pt in range(0,Npts):

        s = S0 + pt * ds        
        s_axis.append(s) 

        for dof in range(0,ndof):
        

            t = (s - s0)/(s1 - s0)

            r = R0.get(dof, 0) + (R1.get(dof, 0) - R0.get(dof, 0)) *  t

            R.set(dof, pt, r)

    # Create the xyz file
    xyz = ""
    for pt in range(0,Npts):
        xyz = xyz + "%i\n" % (nat)
        xyz = xyz + "Frame %i\n" % (pt)
        for at in range(0,nat):
            x = R.get(3*at+0, pt) / units.Angst
            y = R.get(3*at+1, pt) / units.Angst
            z = R.get(3*at+2, pt) / units.Angst
            xyz = xyz + "%s  %8.5f  %8.5f  %8.5f \n" % (E[at], x, y, z)
 
    return R, xyz, s_axis
 