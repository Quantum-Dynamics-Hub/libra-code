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
.. module:: hippynn_methods
   :platform: Unix, Windows
   :synopsis: This module implements the interface to the hippynn package: it exposes a
       trained HIP-NN machine-learning model as a Libra model Hamiltonian, so that NA-MD
       can be driven by a neural network instead of an analytical potential or an
       electronic-structure package.

.. moduleauthor:: Victor M. Freixas

The model is evaluated through the hippynn ``Predictor`` API. The user points this
module at a trained hippynn model and names the output node ("label") that carries
the electronic Hamiltonian; nuclear derivatives are obtained by automatic
differentiation, so no derivative labels need to be trained or stored.

Two kinds of labels are supported, distinguished automatically by the shape of the
predicted quantity at a single geometry:

    * **(nstates, nstates)** - a diabatic Hamiltonian matrix. Off-diagonal elements
      are the diabatic couplings, so non-adiabatic dynamics is fully defined. This
      is the recommended target for NA-MD.
    * **(nstates,)** - a vector of state energies, placed on the diagonal of
      ``ham_dia``. The off-diagonal elements are then zero, i.e. the states are
      uncoupled; useful for adiabatic/Born-Oppenheimer dynamics or as a starting
      point, but it does not by itself describe non-adiabatic transitions.

Requires the ``hippynn`` and ``torch`` packages, which are imported lazily so that
the rest of ``libra_py.models`` is unaffected if they are not installed.

Adapting this module to another hippynn model
---------------------------------------------
Nothing here is specific to a given molecule or model; the molecule enters only
through ``params``. To use your own model you need to set:

    * ``params["model_path"]``   - the directory holding your trained model
    * ``params["label"]``        - the name of the node to read (see above)
    * ``params["Z"]``            - the atomic numbers, in the order your model expects
    * ``params["energy_units"]`` / ``params["length_units"]`` - the units your
      model was trained in; conversion to atomic units happens here.

Two practical notes:

    1. A hippynn model is *not* self-contained: if it was built with custom node
       classes, those classes must be importable whenever the saved graph is
       deserialized (e.g. the defining package must be on ``sys.path``). This applies
       equally to ``params["model_path"]`` and to building the object yourself, since
       both unpickle the same graph. Use ``params["predictor"]`` to pass an
       already-constructed ``GraphModule`` or ``Predictor`` when the model does not
       live in a standard checkpoint directory, or is already in memory - it changes
       *how* the network is obtained, not what its definition requires.
    2. The input node names default to ``"Z"`` (atomic numbers) and ``"R"``
       (positions). If your model names them differently, set ``params["input_Z"]``
       and ``params["input_R"]``.

"""

import os
import sys

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units


class tmp:
    pass


# Building a Predictor is expensive (it loads the network), while Libra calls the
# model function once per trajectory per step. Predictors are therefore cached and
# reused; the key identifies the model, not the geometry.
_PREDICTOR_CACHE = {}


def _get_predictor(params):
    """Return a (predictor, label) pair for the model described by ``params``, cached.

    Args:
        params ( dictionary ): see :func:`hippynn_model`

    Returns:
        tuple: ( hippynn.graphs.Predictor, string ): the predictor and the label to read

    """
    label = params["label"]

    # An already-built Predictor or GraphModule takes precedence over a path: this is
    # the escape hatch for models whose checkpoints cannot be loaded standalone.
    obj = params["predictor"]
    if obj is not None:
        key = id(obj)
        if key not in _PREDICTOR_CACHE:
            from hippynn.graphs import GraphModule, Predictor
            if isinstance(obj, Predictor):
                _PREDICTOR_CACHE[key] = obj      # already built by the user; used as given
            else:
                node = obj.node_from_name(label)
                _PREDICTOR_CACHE[key] = Predictor.from_graph(
                    GraphModule(obj.input_nodes, [node]), requires_grad=True,
                    model_device=params["device"])
        return _PREDICTOR_CACHE[key], label

    model_path = params["model_path"]
    key = (os.path.abspath(model_path), label, params["device"])
    if key not in _PREDICTOR_CACHE:
        from hippynn.graphs import GraphModule, Predictor
        from hippynn.experiment.serialization import load_model_from_cwd

        # Only the model is loaded, never the optimizer/database state that a full
        # checkpoint also carries: dynamics is pure inference. hippynn resolves its
        # files relative to the working directory, hence the chdir.
        cwd = os.getcwd()
        os.chdir(model_path)
        try:
            model = load_model_from_cwd(model_device=params["device"])
        finally:
            os.chdir(cwd)

        # The label may be an interior node of the graph rather than one of its
        # outputs, so it is requested explicitly. Only the sub-graph feeding that node
        # is retained: a trained model also carries the other heads used for fitting
        # (state gradients, couplings, masks, ...), and re-evaluating them at every
        # dynamics step is pure overhead.
        node = model.node_from_name(label)
        _PREDICTOR_CACHE[key] = Predictor.from_graph(
            GraphModule(model.input_nodes, [node]), requires_grad=True,
            model_device=params["device"])

    return _PREDICTOR_CACHE[key], label


def hippynn_model(q, params, full_id):
    """

    A Libra model Hamiltonian backed by a trained hippynn neural network.

    The network is evaluated at the current nuclear geometry; the requested label is
    interpreted either as a diabatic Hamiltonian matrix or as a vector of state
    energies (see the module docstring), and its nuclear derivatives are computed by
    automatic differentiation.

    Args:
        q ( MATRIX(ndof, ntraj) ): coordinates of the particles, in Bohr.
            ``ndof = 3 * natoms``, ordered as x1, y1, z1, x2, y2, z2, ...
        params ( dictionary ): model parameters

            * **params["model_path"]** ( string ): directory with the trained hippynn
                model, as saved by hippynn's experiment machinery. Only the files
                needed to rebuild the network are read (the structure file and the
                model weights); optimizer/database state is not required
                [ required, unless **params["predictor"]** is given ]
            * **params["Z"]** ( list of ints ): atomic numbers of the atoms, in the
                same order as the coordinates in ``q`` [ required ]
            * **params["label"]** ( string ): name of the hippynn node to read
                [ default: "W" ]
            * **params["predictor"]** ( Predictor or GraphModule ): an already-built
                hippynn object to use instead of loading from disk [ default: None ]
            * **params["device"]** ( string ): torch device for the network, e.g.
                "cpu" or "cuda" [ default: "cpu" ]
            * **params["energy_units"]** ( string ): units of the model's output,
                one of "eV", "Ha" [ default: "eV" ]
            * **params["length_units"]** ( string ): length units the model expects,
                one of "Angstrom", "Bohr" [ default: "Angstrom" ]
            * **params["input_Z"]** ( string ): name of the model's atomic-number
                input [ default: "Z" ]
            * **params["input_R"]** ( string ): name of the model's position input
                [ default: "R" ]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(nstates,nstates) ): diabatic Hamiltonian [ units: Ha ]
            * obj.ovlp_dia ( CMATRIX(nstates,nstates) ): overlap of the diabatic states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(nstates,nstates) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinates
                [ units: Ha/Bohr ]
            * obj.dc1_dia ( list of ndof CMATRIX(nstates,nstates) objects ):
                derivative couplings in the diabatic basis [ zero: the diabatic states
                are taken to be strictly diabatic ]

    """

    critical_params = ["Z"]
    default_params = {"model_path": None, "label": "W", "predictor": None,
                      "device": "cpu", "energy_units": "eV", "length_units": "Angstrom",
                      "input_Z": "Z", "input_R": "R"}
    comn.check_input(params, default_params, critical_params)

    if params["model_path"] is None and params["predictor"] is None:
        raise ValueError("Hippynn model: provide either params['model_path'] or params['predictor']")

    import numpy as np
    import torch

    predictor, label = _get_predictor(params)

    Z = np.asarray(params["Z"], dtype=np.int64)
    natoms = Z.shape[0]
    ndof = 3 * natoms

    # unit conversion: the network's units -> atomic units
    e_conv = {"eV": units.ev2Ha, "Ha": 1.0}[params["energy_units"]]
    # positions go the other way: atomic units (Libra) -> the network's units
    r_conv = {"Angstrom": 1.0 / units.Angst, "Bohr": 1.0}[params["length_units"]]

    indx = Cpp2Py(full_id)[-1]
    R = np.array([q.get(i, indx) for i in range(ndof)], dtype=np.float64).reshape(1, natoms, 3)

    z_in = torch.as_tensor(Z[None, :], device=params["device"])
    r_in = torch.tensor(R * r_conv, dtype=torch.float32,
                        device=params["device"], requires_grad=True)

    out = predictor(**{params["input_Z"]: z_in, params["input_R"]: r_in})
    val = out[label] if label in out else out[[k for k in out if getattr(k, "name", None) == label][0]]
    val = val[0]

    # (nstates, nstates) -> full diabatic Hamiltonian; (nstates,) -> energies on the diagonal
    is_matrix = (val.dim() == 2)
    nstates = val.shape[0]

    ham_dia = CMATRIX(nstates, nstates)
    ovlp_dia = CMATRIX(nstates, nstates)
    ovlp_dia.identity()
    d1ham_dia = CMATRIXList()
    dc1_dia = CMATRIXList()
    for i in range(ndof):
        d1ham_dia.append(CMATRIX(nstates, nstates))
        dc1_dia.append(CMATRIX(nstates, nstates))

    # values and their nuclear gradients. Only the unique elements are differentiated;
    # a diabatic Hamiltonian is symmetric, so the lower triangle is filled by copy.
    v = val.detach().cpu().numpy()
    for m in range(nstates):
        for n in range(m, nstates if is_matrix else m + 1):
            elem = val[m, n] if is_matrix else val[m]
            value = (v[m, n] if is_matrix else v[m]) * e_conv
            ham_dia.set(m, n, value * (1.0 + 0.0j))
            if n != m:
                ham_dia.set(n, m, value * (1.0 + 0.0j))

            # d(element)/dR; retain the graph so the remaining elements can be done too
            grad = torch.autograd.grad(elem, r_in, retain_graph=True, allow_unused=True)[0]
            if grad is None:
                continue                       # element is structurally independent of R
            g = grad[0].detach().cpu().numpy().reshape(ndof) * e_conv * r_conv

            for i in range(ndof):
                d1ham_dia[i].set(m, n, g[i] * (1.0 + 0.0j))
                if n != m:
                    d1ham_dia[i].set(n, m, g[i] * (1.0 + 0.0j))

    obj = tmp()
    obj.ham_dia = ham_dia
    obj.ovlp_dia = ovlp_dia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
