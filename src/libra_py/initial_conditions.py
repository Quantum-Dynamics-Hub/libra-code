# *********************************************************************************
# * Copyright (C) 2026 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

import numpy as np

def remove_translation(p):
    """
    Remove total linear momentum.

    Parameters
    ----------
    p : (ndof, ntraj)

    Returns
    -------
    p_clean : (ndof, ntraj)
    """

    ndof, ntraj = p.shape
    natoms = ndof // 3

    P = p.astype(float).reshape(natoms, 3, ntraj).copy()

    for itraj in range(ntraj):

        Ptot = np.sum(P[:, :, itraj], axis=0)

        # distribute equally across atoms
        P[:, :, itraj] -= Ptot / natoms

    return P.reshape(ndof, ntraj)



def remove_rotation(q, p, masses):
    """
    Remove total angular momentum.

    Parameters
    ----------
    q : (ndof, ntraj)
    p : (ndof, ntraj)
    masses : (ndof,)

    Returns
    -------
    p_clean : (ndof, ntraj)
    """

    ndof, ntraj = p.shape
    natoms = ndof // 3

    m = masses.reshape(natoms, 3)[:, 0]

    R = q.astype(float).reshape(natoms, 3, ntraj)
    P = p.astype(float).reshape(natoms, 3, ntraj).copy()

    for itraj in range(ntraj):

        r = R[:, :, itraj].copy()
        mom = P[:, :, itraj].copy()

        # shift coordinates to COM
        r_com = np.sum(m[:, None] * r, axis=0) / np.sum(m)
        r -= r_com

        # total angular momentum
        L = np.sum(np.cross(r, mom), axis=0)

        # inertia tensor
        I = np.zeros((3, 3))

        for i in range(natoms):

            ri = r[i]

            I += m[i] * (
                np.dot(ri, ri) * np.eye(3)
                - np.outer(ri, ri)
            )

        omega = np.linalg.pinv(I) @ L

        # subtract rotational momentum
        for i in range(natoms):

            mom[i] -= m[i] * np.cross(omega, r[i])

        P[:, :, itraj] = mom

    return P.reshape(ndof, ntraj)



def cleanup_momenta(q, p, masses,
                    stages=("translation",
                            "rotation",
                            "translation")):
    """
    Sequential cleanup of momenta.

    Example:
        translation -> rotation -> translation
    """

    p_clean = p.copy()

    for stage in stages:

        if stage == "translation":
            p_clean = remove_translation(p_clean)

        elif stage == "rotation":
            p_clean = remove_rotation(q, p_clean, masses)

        else:
            raise ValueError(f"Unknown stage: {stage}")

    return p_clean


