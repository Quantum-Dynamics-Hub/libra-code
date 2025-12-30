# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for reading and handling xyz files
       This package contains the following functions:
           * cryst2cart(a1,a2,a3,r)

.. moduleauthor:: Alexey V. Akimov

"""

import numpy as np
import libra_py.units as units


def read_trajectory(
    filename,
    scaling = 1.0
):
    """
    Read an MD trajectory from an XYZ file using NumPy.

    Parameters
    ----------
    filename : str
        Path to the XYZ file containing the MD trajectory.
        Standard XYZ format is assumed:
            nat
            comment
            element  x  y  z
            ...

    scaling : float, optional
        How much to scale the raw input data

    Returns
    -------
    R : np.ndarray, shape (3*nat, nsteps)
        Cartesian coordinates of all atoms for all time steps,
        ordered as (x1, y1, z1, x2, y2, z2, ...).

    E : list of str, length nat
        Atomic symbols (element names), in the order they appear
        in the XYZ file.

    Notes
    -----
    - All frames in the XYZ file are read.
    - The array layout (3*nat, nsteps) matches legacy Libra-style
      coordinate storage and is convenient for mass-weighting
      and correlation-matrix analysis.

    Examples
    --------
    Read an XYZ trajectory in Angstrom and return coordinates in Angstrom
    (default behavior):

    >>> R, E = read_trajectory("traj.xyz")

    Read an XYZ trajectory in Angstrom and convert them to Bohr:

    >>> R, E = read_trajectory(
    ...     "traj.xyz",
    ...     scaling = units.Angst
    ... )


    Convert the returned array to (nsteps, nat, 3) form:

    >>> nat = len(E)
    >>> nsteps = R.shape[1]
    >>> R3 = R.T.reshape(nsteps, nat, 3)

    Access the trajectory of atom i (0-based indexing):

    >>> i = 0
    >>> ri = R[3*i:3*i+3, :]   # shape (3, nsteps)
    """

    # -----------------------------
    # Read file
    # -----------------------------
    with open(filename, "r") as f:
        lines = f.readlines()

    nat = int(lines[0].split()[0])
    block_size = nat + 2
    nsteps = len(lines) // block_size

    # -----------------------------
    # Allocate arrays
    # -----------------------------
    R = np.zeros((3 * nat, nsteps), dtype=float)
    E = []

    # -----------------------------
    # Parse trajectory
    # -----------------------------
    for t in range(nsteps):
        offset = t * block_size + 2

        for i in range(nat):
            tokens = lines[offset + i].split()
            name = tokens[0]
            xyz = np.array(tokens[1:4], dtype=float) * scaling

            R[3*i:3*i+3, t] = xyz

            if t == 0:
                E.append(name)

    return R, E



def names_to_masses(E, PT, dtype=float):
    """
    Construct a Cartesian mass vector of length 3N from atomic element names.

    Parameters
    ----------
    E : sequence of str, length N
        Atomic symbols (element names), e.g. ["Ti", "O", "O", ...].

    PT : dict[str, float]
        Mapping from element name to atomic mass.
        Masses are assumed to be given in atomic mass units (amu).

    dtype : data-type, optional
        Data type of the returned array. Default is float.

    Returns
    -------
    M : np.ndarray, shape (3*N,)
        Cartesian mass vector, where each atomic mass is repeated
        three times:
            (m1, m1, m1, m2, m2, m2, ..., mN, mN, mN)

    Raises
    ------
    KeyError
        If an element in `E` is not found in the mass dictionary `PT`.

    Examples
    --------
    >>> PT = {"Ti": 47.9, "O": 16.0}
    >>> E = ["Ti", "O", "O"]
    >>> M = names_to_masses(E, PT)
    >>> M
    array([47.9, 47.9, 47.9, 16. , 16. , 16. , 16. , 16. , 16. ])

    Use the mass vector for mass-weighting Cartesian coordinates:

    >>> Rmw = R * np.sqrt(M)[:, None]
    """

    try:
        masses = np.array([PT[e] for e in E], dtype=dtype)
    except KeyError as exc:
        raise KeyError(f"Mass for element '{exc.args[0]}' not found in PT") from None

    # Repeat each atomic mass three times (x, y, z)
    M = np.repeat(masses, 3)

    return M



def write_trajectory(filename, R, E, scaling = 1.0):
    """
    Write a trajectory to an XYZ-format file.

    Args:
        filename : str
            Path to the output XYZ file.
        R : np.ndarray of shape (3*N, nsteps)
            Cartesian coordinates for all atoms and timesteps.
        E : list of str, length N
            Atom labels (elements) for each atom.
        scaling : float, default 1.0
            How much scale the input data before writing to file

    Usage:
        >>> R, E = read_trajectory("coords.xyz", scaling = units.Angst) # Angstrom to Bohr
        >>> write_trajectory("out.xyz", R, E, scaling = 1.0/units.Angst) # Bohr to Angstrom
    """
    R = np.array(R)
    N = len(E)
    nsteps = R.shape[1]

    # -----------------------------
    # Write XYZ trajectory
    # -----------------------------
    with open(filename, 'w') as f:
        for t in range(nsteps):
            f.write(f"{N}\n")
            f.write(f"Step {t+1}\n")
            for i in range(N):
                x = R[3*i + 0, t] * scaling
                y = R[3*i + 1, t] * scaling
                z = R[3*i + 2, t] * scaling
                f.write(f"{E[i]:<3} {x:>15.8f} {y:>15.8f} {z:>15.8f}\n")




