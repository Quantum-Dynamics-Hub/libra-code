# *********************************************************************************
# * Copyright (C) 2024 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: md_align
   :platform: Unix, Windows
   :synopsis: this module implements functions for removing COM motion and overall
       rotation of the system

       Note: This trajectory alignment code is a Python translation of the Fortran code of Sean Fisher found here:
       https://socks.lbl.gov/jdbachan/nwchemxx/-/blob/5da5730af2a818947db1f9c621b36315f59a8eed/contrib/qmd_tools/qmd_analysis.f90

       Example:

       >>> align_trajectory("C20-md.xyz", "C20-md-aligned.xyz")

.. moduleauthor:: Alexey V. Akimov

"""


import numpy as np

from liblibra_core import *
import libra_py.packages.qe.methods as qe
from . import data_conv
from . import units


def compute_com(R, M):
    """
    Args:

        R (MATRIX(3N, 1)): coordinates of N atoms
        M (list of N doubles)): masses of N atoms

    Returns:
        tuple: (X, Y, Z): coordinates of the COM

    """
    totM = sum(M)
    ndof = R.num_of_rows
    nat = int(ndof / 3)
    X, Y, Z = 0.0, 0.0, 0.0

    for i in range(nat):
        X = X + (M[i] / totM) * R.get(3 * i + 0, 0)
        Y = Y + (M[i] / totM) * R.get(3 * i + 1, 0)
        Z = Z + (M[i] / totM) * R.get(3 * i + 2, 0)

    return X, Y, Z


def remove_com(R, M):
    """
    Args:

        R (MATRIX(3N, 1)): coordinates of N atoms
        M (list of N doubles)): masses of N atoms

    Returns:
        MATRIX(3N,1): molecular coordinates with subtracted COM

    """

    X, Y, Z = compute_com(R, M)
    R_centered = MATRIX(R)

    ndof = R.num_of_rows
    nat = int(ndof / 3)

    for i in range(nat):
        R_centered.add(3 * i + 0, 0, -X)
        R_centered.add(3 * i + 1, 0, -Y)
        R_centered.add(3 * i + 2, 0, -Z)
    return R_centered


def compute_covar(coord, ref):
    """
    Args:
        coord (MATRIX(3N, 1)): coordinates of N atoms
        ref  (MATRIX(3N, 1)): reference structure of N atoms

    Returns:
        MATRIX(3,3): covariance matrix
    """

    R = MATRIX(3, 3)
    ndof = ref.num_of_rows
    nat = int(ndof / 3)

    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            res = 0.0
            for k in range(nat):
                res = res + coord.get(3 * k + i, 0) * ref.get(3 * k + j, 0)
            R.set(i, j, res)
    return R


def compute_rotation_mtx(R):
    """
    Args:
        R (MATRIX(3,3)): covariance matrix

    Returns:
        MATRIX(3,3): rotation matrix
    """

    F = MATRIX(4, 4)
    F.set(0, 0, R.get(0, 0) + R.get(1, 1) + R.get(2, 2))
    F.set(1, 0, R.get(1, 2) - R.get(2, 1))
    F.set(0, 1, F.get(1, 0))
    F.set(2, 0, R.get(2, 0) - R.get(0, 2))
    F.set(0, 2, F.get(2, 0))
    F.set(3, 0, R.get(0, 1) - R.get(1, 0))
    F.set(0, 3, F.get(3, 0))
    F.set(1, 1, R.get(0, 0) - R.get(1, 1) - R.get(2, 2))
    F.set(2, 1, R.get(0, 1) + R.get(1, 0))
    F.set(1, 2, F.get(2, 1))
    F.set(3, 1, R.get(0, 2) + R.get(2, 0))
    F.set(1, 3, F.get(3, 1))
    F.set(2, 2, R.get(1, 1) - R.get(0, 0) - R.get(2, 2))
    F.set(3, 2, R.get(1, 2) + R.get(2, 1))
    F.set(2, 3, F.get(3, 2))
    F.set(3, 3, R.get(2, 2) - R.get(0, 0) - R.get(1, 1))

    #
    # INTEGER :: LDA=4
    # INTEGER :: LWORK=3*4-1
    # INTEGER :: INFO
    # CHARACTER, PARAMETER :: JOBZ = 'V'  # Compute eigenvalues and eigenvectors
    # CHARACTER, PARAMETER :: UPLO = 'U'  # Upper triangle of A is stored;
    # DOUBLE PRECISION,DIMENSION(3*4-1) :: WORK
    # call DSYEV(JOBZ,UPLO,4,F,LDA,val,WORK,LWORK,INFO)

    # solve_eigen(MATRIX& H, MATRIX& E, MATRIX& C, int symm);
    E = MATRIX(4, 4)
    C = MATRIX(4, 4)
    solve_eigen(F, E, C, 0)

    tmp = data_conv.MATRIX2nparray(C.col(3), np.float64)
    q = [tmp[0, 0], tmp[1, 0], tmp[2, 0], tmp[3, 0]]
    # print(q)

    U = MATRIX(3, 3)
    U.set(0, 0, q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2)
    U.set(1, 1, q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2)
    U.set(2, 2, q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2)

    U.set(1, 0, 2.0 * (q[1] * q[2] + q[0] * q[3]))
    U.set(0, 1, 2.0 * (q[1] * q[2] - q[0] * q[3]))
    U.set(2, 0, 2.0 * (q[1] * q[3] - q[0] * q[2]))
    U.set(0, 2, 2.0 * (q[1] * q[3] + q[0] * q[2]))
    U.set(2, 1, 2.0 * (q[2] * q[3] + q[0] * q[1]))
    U.set(1, 2, 2.0 * (q[2] * q[3] - q[0] * q[1]))

    return U


def remove_rotation(R, U):
    """
    Args:

        R (MATRIX(3N, 1)): coordinates of N atoms
        U (MATRIX(3,3))): roation of the current coordinate to the reference

    Returns:
        MATRIX(3N,1): molecular coordinates with removed rotation

    """

    ndof = R.num_of_rows
    nat = int(ndof / 3)

    r = MATRIX(3, 1)

    R_rotated = MATRIX(R)

    for i in range(nat):
        r.set(0, 0, R.get(3 * i + 0, 0))
        r.set(1, 0, R.get(3 * i + 1, 0))
        r.set(2, 0, R.get(3 * i + 2, 0))

        r = U * r

        R_rotated.set(3 * i + 0, 0, r.get(0, 0))
        R_rotated.set(3 * i + 1, 0, r.get(1, 0))
        R_rotated.set(3 * i + 2, 0, r.get(2, 0))

    return R_rotated


def align_trajectory(in_filename, out_filename):
    # Read in the original trajectory
    R, E = qe.read_md_data_xyz2(in_filename, {})
    # Atomic weights of elements from atomic number 1 to 80
    atomic_weights = {
        "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
        "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
        "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
        "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
        "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
        "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
        "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.971, "Br": 79.904,
        "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
        "Nb": 92.906, "Mo": 95.95, "Tc": 98, "Ru": 101.07, "Rh": 102.91,
        "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
        "Sb": 121.76, "Te": 127.6, "I": 126.9, "Xe": 131.29, "Cs": 132.91,
        "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
        "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
        "Dy": 162.5, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05,
        "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21,
        "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
        "Tl": 204.38, "Pb": 207.2, "Bi": 208.98
    }

    # print(E)
    # By default, the reading converts the coordinates to Bohrs, but
    # here, we just want to keep the original Angstrom units, so:
    R.scale(-1, -1, 1.0 / units.Angst)

    # Masses
    M = [atomic_weights[e] for e in E]

    # print( R.num_of_rows, R.num_of_cols, len(M) )
    ndof = R.num_of_rows
    nat = int(ndof / 3)
    nsteps = R.num_of_cols

    f = open(out_filename, "w")
    R0 = R.col(0)
    ref = remove_com(R0, M)
    for i in range(nsteps):
        # print(i)

        rr = None
        if i == 0:
            rr = ref
        else:
            r = R.col(i)
            r = remove_com(r, M)
            covar = compute_covar(r, ref)
            U = compute_rotation_mtx(covar)
            r = remove_rotation(r, U)
            rr = r

        res = "%5i\nstep%5i\n" % (nat, i)

        for n in range(nat):
            x, y, z = rr.get(3 * n, 0), rr.get(3 * n + 1, 0), rr.get(3 * n + 2, 0)

            res = res + "%s  %5.3f  %5.3f  %5.3f\n" % (E[n], x, y, z)
        f.write(res)
    f.close()
