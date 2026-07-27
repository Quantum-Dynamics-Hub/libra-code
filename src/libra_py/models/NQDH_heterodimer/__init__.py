# *********************************************************************************
# * Copyright (C) 2026 Victor M. Freixas
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: NQDH_heterodimer
   :platform: Unix, Windows
   :synopsis: A trained neural-network Hamiltonian for a halofluorescein AB heterodimer.
       Like the LVC or Shin-Metiu models, this describes one specific system, but its
       parameters are the weights of a machine-learning model rather than a handful of
       constants.

.. moduleauthor:: Victor M. Freixas

The system is a 122-atom halofluorescein dimer: two chromophores whose lowest excited
states are delocalized combinations of excitations on either half. The model predicts
the diabatic electronic Hamiltonian ``W(R)`` of the three lowest singlet states
(S0, S1, S2). Its eigenvalues are the adiabatic energies, and its off-diagonal S1-S2
element is the diabatic coupling that drives transitions between the two excited states.

The ground state is not coupled to the excited states in the reference data, and the
network was built with that structure imposed, so ``W[0, 1] = W[0, 2] = 0`` exactly
rather than approximately.

The model was trained on AM1/CIS reference data (energies, gradients and non-adiabatic
coupling vectors) computed with NEXMD. Propagated with Ehrenfest dynamics over a thermal
ensemble, it reproduces the S1/S2 vibronic beating of the dimer: recurrences with a
period of about 18 fs, decaying over roughly 90 fs.

Usage
-----
The model is evaluated through the hippynn interface in
:mod:`libra_py.packages.hippynn.methods`, which this module wraps with the paths and the
atomic composition already filled in::

    from libra_py.models import NQDH_heterodimer

    params = {"model0": 0, "nstates": 3}
    obj = NQDH_heterodimer.compute_model(q, params, full_id)

``q`` holds the nuclear coordinates in Bohr, ordered x1, y1, z1, x2, y2, z2, ..., for
the 122 atoms in the order given by :func:`get_atomic_numbers`. A set of thermal initial
conditions for the dynamics is provided by :func:`get_initial_conditions`.

Requires ``hippynn`` and ``torch``, which are imported only when the model is evaluated.

Notes
-----
The ``nqdh_nodes`` sub-package holds the custom hippynn node classes the network was
built with. A saved hippynn model records the module path of such classes, so they have
to be importable for the model to load; they are shipped here for that reason and are
not otherwise used. The loss-function nodes among them are included because hippynn
graph nodes are doubly linked, so the saved sub-graph retains references to the nodes
that were attached to it during training.

"""

import os

_HERE = os.path.dirname(os.path.abspath(__file__))

#: the trained network: a hippynn graph mapping (Z, R) to the diabatic Hamiltonian W
GRAPH_FILE = os.path.join(_HERE, "w_graph.pt")

#: thermal initial conditions (positions, velocities, masses) for the dynamics
ICS_FILE = os.path.join(_HERE, "ics.npz")

#: the label of the node holding the diabatic Hamiltonian
LABEL = "W"

_CACHE = {}


def get_atomic_numbers():
    """The atomic numbers of the 122 atoms, in the order the model expects them.

    Returns:
        list of ints: atomic numbers

    """
    import numpy as np
    return np.load(ICS_FILE)["Z"].astype(int).tolist()


def get_initial_conditions():
    """Thermal initial conditions sampled from a ground-state trajectory at 300 K.

    Returns:
        tuple: ( R, V, masses ), where

            * R ( numpy array (ntraj, natoms, 3) ): positions [ units: Bohr ]
            * V ( numpy array (ntraj, natoms, 3) ): velocities [ units: a.u. ]
            * masses ( numpy array (natoms) ): atomic masses [ units: a.u. ]

    """
    import numpy as np
    data = np.load(ICS_FILE)
    return (np.asarray(data["R_bohr"], dtype=float),
            np.asarray(data["V_au"], dtype=float),
            np.asarray(data["masses_au"], dtype=float))


def compute_model(q, params, full_id):
    """

    The diabatic Hamiltonian of the halofluorescein heterodimer, from a trained network.

    This is a thin wrapper around :func:`libra_py.packages.hippynn.methods.hippynn_model`
    with the model file, the label and the atomic composition of this particular system
    already set, so that only the dynamics parameters have to be supplied.

    Args:
        q ( MATRIX(ndof, ntraj) ): coordinates of the atoms, ndof = 3 * 122,
            ordered x1, y1, z1, x2, y2, z2, ... [ units: Bohr ]
        params ( dictionary ): model parameters. Any key accepted by
            :func:`libra_py.packages.hippynn.methods.hippynn_model` may be given to
            override the defaults set here, e.g. ``params["device"] = "cuda"``.
        full_id ( intList ): the trajectory identifier

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(3,3) ): diabatic Hamiltonian [ units: Ha ]
            * obj.ovlp_dia ( CMATRIX(3,3) ): overlap of the diabatic states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(3,3) ): its nuclear derivatives
                [ units: Ha/Bohr ]
            * obj.dc1_dia ( list of ndof CMATRIX(3,3) ): derivative couplings in the
                diabatic basis [ zero ]

    """
    from libra_py.packages.hippynn.methods import hippynn_model

    if "predictor" not in _CACHE:
        import torch
        _CACHE["predictor"] = torch.load(GRAPH_FILE, map_location="cpu", weights_only=False)

    params.setdefault("Z", get_atomic_numbers())
    params.setdefault("label", LABEL)
    params["predictor"] = _CACHE["predictor"]

    return hippynn_model(q, params, full_id)
