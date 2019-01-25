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

import os
import math
import sys
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import units


"""
  ========= Theory ==============
  
  Convention on the normal modes transformations:

  U^T * H * U = W^2       (1)

  x = U * q               (2a)

  so q = U^T * x          (2b)
                       
  x = M^(1/2) * (R-R_ave) (3)

  Here:
  H - mass-scaled Hessian (dynamical matrix): H_ij = (d^2E/dx_i * dx_j )/sqrt(m_i * m_j)
  U, q - normal modes (q are columns of U)
  x - mass-scaled and shifted coordinates  
  R - Cartesian coordinates
      
  Note: to relate to the definitions here (as of 11/7/2018), use: U = D^T
  https://www.pci.uni-heidelberg.de//tc/usr/mctdh/doc/vcham/latex/nmode_coo.pdf

  Also useful:
  https://nanohub.org/courses/FATM/01a/asset/256

"""


def covariance_matrix(X, M, flag):
    """
    Computes the covariance matrix K^x = <sqrt(m_i * m_j) * x_i * x_j }>

    Input:
    X (MATRIX(ndof, nsteps)) - data collected along a trajectory: can be R, V, or A
    M (MATRIX(ndof, 1)) - masses of all DOFs 
    flag (int) - controls how to compute variance
          0 - using the data as they are (no centering)
          1 - using the fluctuations of the data around the mean (do centering)

    Returns:
    K (MATRIX(ndof, ndof)) - the matrix of covariance of all DOFs averaged over the trajectory
    """

    ndof = X.num_of_rows
    nsteps = X.num_of_cols

    # To speed-up the calculations, minimize the number of sqrt calls
    sM = MATRIX(ndof, 1)
    for i in xrange(ndof):
        sM.set(i, math.sqrt(M.get(i)) )

    K = None
    if flag==0:
        K = covariance(X)
    elif flag==1:
        dX = deviation(X)
        K = covariance(dX)

    # Compute the covariance matrix
    for i in xrange(ndof):
        for j in xrange(ndof):
            mij = sM.get(i) * sM.get(j)
            K.scale(i,j, 0.5*mij)

    return K


def visualize_modes(E, R, U, M, w, params):
    """
    This function "visualizes" a particular collective modes for a particular set of data.
    The visualization means we generate an "xyz" file, showing the trajectory with defined number
    of repetitions of the mode.

    Input:
    E (list of ndof/3) - atom names (elements) of all atoms
    R ( MATRIX(ndof x nsteps-1) ) - coordinates of all DOFs for all mid-timesteps
    U ( MATRIX(ndof, ndof) ) - eigenvectors defining the collective modes in terms of the original DOFs, see the theory
    M ( MATRIX(ndof x 1) ) - masses of all DOFs
    w (list of ndof) - frequencies of all modes
    params (dictionary) - parameters controlling how to do the visualization
      Contains keyword-value pairs:
      "scale" (double) - mode amplification factor (for better visualization)
      "print_modes" (list of integers) - indices of the modes to hande. Indexing starts with 0 and should 
      be consistent with other data arrays provided - e.g. "w"
      "prefix" (string) - the name of the prefix of the files, to where the modes are printed out
      "nperiods" (integer) - the number of periods of motion to repeat
      "nsteps" (integer) - how many steps should be in the "visualization" trajectory. Controlls the resolution of the modes

    All quantities are in atomic units

    Returns:
    Simply prints out a number of files with the trajectories (using Angstrom units) showing the
    periodic motion along the modes of interest
 
    """

    ndof = R.num_of_rows
    nat = ndof/3

    scl = params["scale"]

    for mode in params["print_modes"]:
    
        filename = params["prefix"]+"_mode%i.xyz" % (mode)

        dt = 0.0
        if w.get(mode)>0.0:
            dt = 2.0*math.pi * params["nperiods"] / ( w.get(mode) * params["nsteps"])

        f = open(filename, "w") 
        for t in xrange(params["nsteps"]):
            f.write("%i\n" % (nat))
            f.write("step %i = \n" % (t) )
 
            for at in xrange(nat): 
                x = (R.get(3*at+0, 0) + (U.get(3*at+0, mode)/math.sqrt(M.get(3*at+0))) * math.sin(w.get(mode) * t*dt) * scl)/units.Angst
                y = (R.get(3*at+1, 0) + (U.get(3*at+1, mode)/math.sqrt(M.get(3*at+1))) * math.sin(w.get(mode) * t*dt) * scl)/units.Angst
                z = (R.get(3*at+2, 0) + (U.get(3*at+2, mode)/math.sqrt(M.get(3*at+2))) * math.sin(w.get(mode) * t*dt) * scl)/units.Angst

                f.write("%s  %8.5f %8.5f %8.5f \n" % (E[at], x, y, z) )
        
        f.close()
    
   

def compute_cov(R, V, A, M, E, params):
    """
    Computes and visualizes (as the trajectories) normal modes following the
    methods described in:

    (1) Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics 
    or Monte Carlo Simulation. J. Chem. Phys. 2003, 120, 1-4.

    Input:
    R ( MATRIX(ndof x nsteps-1) ) - coordinates of all DOFs for all mid-timesteps
    V ( MATRIX(ndof x nsteps-1) ) - velocities of all DOFs for all mid-timesteps
    A ( MATRIX(ndof x nsteps-1) ) - accelerations of all DOFs for all mid-timesteps
    M ( MATRIX(ndof x 1) ) - masses of all DOFs
    E (list of ndof/3) - atom names (elements) of all atoms

    params (dictionary) - parameters controlling the computations, including the visualization
    (see the visualize_modes(E, R, U, w, params) description)

    In addition:
    params["verbosity"]  [int]    - level to control verbosity
    params["visualize"]  [int]    - flag to control whether we want to produce additional files (with normal modes)
                                   0 - not to
                                   1 - do it

    All quantities are in atomic units
 
    """
    verbosity = params["verbosity"]

    if verbosity>0:
        print "========= Normal modes calculations according to: ============================="    
        print "Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics\
        or Monte Carlo Simulation. J. Chem. Phys. 2003, 120, 1-4.\n"

    ndof = R.num_of_rows
    nat = ndof/3
    cov_flag = params["cov_flag"]
    
    if verbosity>0:
        print "Computing covariance matrix of positions\n"; 
    K_r = covariance_matrix(R, M, cov_flag)
    if verbosity>0:
        print "Computing covariance matrix of velocities\n"; 
    K_v = covariance_matrix(V, M, cov_flag)
    if verbosity>0:
        print "Computing covariance matrix of accelerations\n"; 
    K_a = covariance_matrix(A, M, cov_flag)

    w_r = MATRIX(ndof, ndof);  U_r = MATRIX(ndof, ndof)
    w_v = MATRIX(ndof, ndof);  U_v = MATRIX(ndof, ndof)
    w_a = MATRIX(ndof, ndof);  U_a = MATRIX(ndof, ndof)

    if verbosity>0:
        print "Eigenvalue solver for covariance matrix of positions\n"
    solve_eigen(K_r, w_r, U_r, 0)
    if verbosity>0:
        print "Eigenvalue solver for covariance matrix of velocities\n"
    solve_eigen(K_v, w_v, U_v, 0)
    if verbosity>0:
        print "Eigenvalue solver for covariance matrix of accelerations\n"
    solve_eigen(K_a, w_a, U_a, 0)

    if verbosity>1:
        print "K_r:"; K_r.show_matrix()
        print "K_r eigenvalues:";  w_r.show_matrix()
        print "K_r eigenvectors:"; U_r.show_matrix()

        print "K_v:"; K_v.show_matrix()
        print "K_v eigenvalues:";  w_v.show_matrix()
        print "K_v eigenvectors:"; U_v.show_matrix()

        print "K_a:"; K_a.show_matrix()
        print "K_a eigenvalues:";  w_a.show_matrix()
        print "K_a eigenvectors:"; U_a.show_matrix()


    w = MATRIX(ndof, 1)      
    for dof in xrange(ndof):
        if w_r.get(dof, dof)>0.0:
            w.set(dof, 0,  math.sqrt( math.fabs( w_v.get(dof, dof)/w_r.get(dof, dof)) )  )
        else:
            w.set(dof, 0, 0.0)
    w_inv_cm = w / units.inv_cm2Ha
    if verbosity>0:
        print "Angular frequencies (derived from w_v/w_r)"
        w_inv_cm.show_matrix()


    w2 = MATRIX(ndof, 1)      
    for dof in xrange(ndof):
        if w_r.get(dof, dof)>0.0:
            w2.set(dof, 0,  math.pow( math.fabs( w_a.get(dof, dof)/w_r.get(dof, dof) ) , 0.25)  )
        else:
            w2.set(dof, 0, 0.0)
    w2_inv_cm = w2 / units.inv_cm2Ha
    if verbosity>0:
        print "Angular frequencies (derived from w_a/w_r)"
        w2_inv_cm.show_matrix()

    if params["visualize"]==1:
        if verbosity>0:
            print "Visualizing modes based on velocities covariance\n" 
        prefix = params["prefix"]
        params.update({"prefix": prefix+"_velocity"})
        visualize_modes(E, R, U_v, M, w, params);

        if verbosity>0:
            print "Visualizing modes based on accelerations covariance\n" 
        params.update({"prefix": prefix+"_acceleration"})
        visualize_modes(E, R, U_a, M, w2, params);

    if verbosity>0:
        print "========= Done with the Normal modes calculations =============================" 

    return w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a


def compute_cov2(R, A, M, E, T, params):
    """
    Computes and visualizes (as the trajectories) normal modes following the
    methods described in:

    (1) Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix:
    Application to Crystalline Explosives PETN and RDX. J. Chem. Phys. 2015, 142, 134110.

    Input:
    R ( MATRIX(ndof x nsteps-1) ) - coordinates of all DOFs for all mid-timesteps
    A ( MATRIX(ndof x nsteps-1) ) - accelerations of all DOFs for all mid-timesteps
    M ( MATRIX(ndof x 1) ) - masses of all DOFs
    E (list of ndof/3) - atom names (elements) of all atoms
    T (double) - temperature of simulation (in K)

    params (dictionary) - parameters controlling the computations, including the visualization
    (see the visualize_modes(E, R, U, w, params) description)

    In addition:
    params["verbosity"]  [int]    - level to control verbosity
    params["visualize"]  [int]    - flag to control whether we want to produce additional files (with normal modes)
                                   0 - not to
                                   1 - do it

    All quantities are in atomic units
 
    """

    verbosity = params["verbosity"]

    if verbosity>0:
        print "========= Normal modes calculations according to: ============================="
        print "Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix:\
        Application to Crystalline Explosives PETN and RDX. J. Chem. Phys. 2015, 142, 134110.\n" 
    
    ndof = R.num_of_rows
    nat = ndof/3
    cov_flag = params["cov_flag"]

    if verbosity>0:
        print "Computing covariance matrix of accelerations\n"
    K_a = None
    if cov_flag==0:
        K_a = covariance(A)
    elif cov_flag==1:
        dA = deviation(A)
        K_a = covariance(dA)

    k = units.boltzmann / units.hartree
    K_a *= (1.0/(k*T))

    if verbosity>0:
        print "Eigenvalue solver for covariance matrix of accelerations\n";
    w_a = MATRIX(ndof, ndof);  U_a = MATRIX(ndof, ndof)
    solve_eigen(K_a, w_a, U_a, 0)

    if verbosity>1:
        print "K_a:"; K_a.show_matrix()
        print "K_a eigenvalues:";  w_a.show_matrix()
        print "K_a eigenvectors:"; U_a.show_matrix()

    w = MATRIX(ndof, 1)      
    for dof in xrange(ndof):
        if w_a.get(dof, dof)>0.0:
            w.set(dof, 0,  math.sqrt(w_a.get(dof, dof)) ) 
    w_inv_cm = w / units.inv_cm2Ha
    if verbosity>0:
        print "Frequencies (cm^-1)";  w_inv_cm.show_matrix()

    if params["visualize"]>0:
        if verbosity>0:
            print "Visualizing modes based on accelerations covariance matrix\n" 
        visualize_modes(E, R, U_a, M, w, params);

    if verbosity>0:
        print "========= Done with the Normal modes calculations ============================="

    return w_a, w_inv_cm, U_a


def compute_dynmat(R, D, M, E, params):
    """
    Computes and visualizes (as the trajectories) normal modes 
    using the dynamic matrix:   D_ij = [1/sqrt(m_i * m_j)]  d^2E/dR_i dR_j

    Here, H_ij = d^2E/dR_i dR_j is the Hessian

    Input:
    R ( MATRIX(ndof x nsteps-1) ) - coordinates of all DOFs for all mid-timesteps
    S ( MATRIX(ndof x nsteps-1) ) - the dynamic matrix
    E (list of ndof/3) - atom names (elements) of all atoms

    params (dictionary) - parameters controlling the computations, including the visualization
    (see the visualize_modes(E, R, U, w, params) description)

    In addition:
    params["verbosity"]  [int]    - level to control verbosity
    params["visualize"]  [int]    - flag to control whether we want to produce additional files (with normal modes)
                                   0 - not to
                                   1 - do it

    All quantities are in atomic units
 
    """

    verbosity = params["verbosity"]

    if verbosity>0:
        print "========= Normal modes calculations using the provided dynamical matrix ============================="
    
    ndof = M.num_of_rows
    nat = ndof/3

    w_a = MATRIX(ndof, ndof);  U_a = MATRIX(ndof, ndof)

    solve_eigen(D, w_a, U_a, 0)

    if verbosity>1:
        print "Dynamic matrix:"; D.show_matrix()
        print "Its eigenvalues:";  w_a.show_matrix()
        print "Its eigenvectors:"; U_a.show_matrix()


    w = MATRIX(ndof, 1)      
    for dof in xrange(ndof):
        if w_a.get(dof, dof)>0.0:
            w.set(dof, 0,  math.sqrt(w_a.get(dof, dof)) )            
    w_inv_cm = w / units.inv_cm2Ha
    if verbosity>0:
        print "Frequencies (cm^-1)";  w_inv_cm.show_matrix()

    if params["visualize"]>0:
        if verbosity>0:
            print "Visualizing modes based on dynamic matrix\n" 
        visualize_modes(E, R, U_a, M, w, params);

    if verbosity>0:
        print "========= Done with the Normal modes calculations ============================="

    return w, w_inv_cm, U_a
 


def get_xyz(E, R, M, U, mode):
    """
    This function returns a string in the xyz format with X, Y, Z and UX, UY, UZ
    where X,Y,Z are the coordinates, UX, UY, UZ - vectors coming from those coordinates - e.g. normal modes

    E [list of ndof/3]            - atom names (elements) of all atoms
    R [MATRIX(ndof x nsteps-1) ]  - coordinates of all DOFs for all mid-timesteps
    M [MATRIX(ndof x 1) ]         - masses of all DOFs
    U [MATRIX(ndof x ndof) ]      - a matrix containing normal mode vectors 
    mode [int]                    - index of the normal mode that we want to visualize

    Returns: 
    a string representing an xyz file

    """

    natoms = len(E)
    res = "%3i\nComment\n" % (natoms)
    
    for i in xrange(natoms):
        x,y,z = R.get(3*i,0), R.get(3*i+1,0), R.get(3*i+2,0)        
        ux = U.get(3*i+0, mode)/math.sqrt(M.get(3*i+0))
        uy = U.get(3*i+1, mode)/math.sqrt(M.get(3*i+1))
        uz = U.get(3*i+2, mode)/math.sqrt(M.get(3*i+2))
        
        res = res + "%s  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n" % (E[i], x,y,z, ux,uy,uz)

    return res
