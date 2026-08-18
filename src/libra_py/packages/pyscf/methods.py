# *********************************************************************************
# * Copyright (C) 2026  Jieyang Gu and Alexey V. Akimov
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
   :synopsis: This module implements an adapter between the ABC interface and the compute_adi 

.. moduleauthor::
       Alexey V. Akimov, Jieyang Gu

"""
import numpy as np

from liblibra_core import CMATRIX, CMATRIXList, Cpp2Py

from libra_py import data_conv
from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
)

def _q_to_geometry(q, itraj, atom_labels):
    """Libra stores nuclear coordinates in Bohr; keep them as ``coords_bohr``."""
    coords = q.col(itraj)
    coordinates = data_conv.MATRIX2nparray(coords, np.float64).reshape(-1, 3)
    return MolecularGeometry( atom_labels=tuple(atom_labels), coords_bohr=coordinates )

class tmp:
    pass

# =============================================================================
# Input conversion helpers
# =============================================================================

def _get_trajectory_index(full_id):
    if isinstance(full_id, (list, tuple)):
        return int(full_id[-1])
    Id = Cpp2Py(full_id)
    return int(Id[-1])

def _resolve_gradient_state(params, itraj):
    """Translate ``params["gradient_state"]`` into what ``ES_Request`` expects.

    Accepted values:
        None    -> no gradient
        int     -> gradient for that one state
        "all"   -> gradients for every state
        "active"-> gradient for the trajectory's current active state, read
                   from ``params["act_state"][itraj]`` (written every step by
                   the dynamics driver, e.g. libra_py.dynamics.tsh.compute).

    ``act_state`` is absent on the warm-up calls that precede the propagation
    loop (libra_py.dynamics.tsh.compute calls the model at lines 1018-1019 but
    only fills act_state at line 1082), so a missing entry falls back to the
    ground state rather than raising -- the same convention the mopac, dftbplus,
    and cp2k interfaces use.
    """
    gradient_state = params.get("gradient_state", "all")

    if isinstance(gradient_state, str):
        gradient_state = gradient_state.lower()

    if gradient_state in (None, "all"):
        return gradient_state

    if gradient_state == "active":
        return params.setdefault("act_state", {}).get(itraj, 0)

    return int(gradient_state)

def _params_to_request(params, itraj): #possible fields of params: "nstates", "H_soc",
                                       #gradient_state", "hes", "nacv", "time_overlap"
    nstates = params.get("nstates", 2)

    return ES_Request(
        n_singlets=nstates,
        n_triplets=0,
        H_soc=params.get("H_soc", False),
        gradient_state=_resolve_gradient_state(params, itraj),
        hessian_state=params.get("hessian_state"),
        nacv=params.get("nacv", False),
        time_overlap=params.get("time_overlap", True),
    )

# =============================================================================
# ES_Result -> Libra type conversion helpers
# =============================================================================

def _energies_to_ham_adi(energies):
    return data_conv.nparray2CMATRIX( np.diag(energies) )

def _gradients_to_d1ham_adi(gradients, nstates, natoms):
    """Pack ES_Result.gradients into Libra's ``d1ham_adi``.

    ``gradients`` is one ``(natoms, 3)`` block per electronic state, or None for
    states whose gradient was not requested. Libra wants one
    CMATRIX(nstates, nstates) per nuclear dof, with ``dof = 3 * atom + xyz``,
    carrying dE_state/dq on its diagonal.
    """
    per_dof = np.zeros((3 * natoms, nstates, nstates), dtype=np.complex128)

    for state, gradient in enumerate(gradients):

        if gradient is None:
            continue

        per_dof[:, state, state] = np.asarray(gradient).reshape(-1)

    d1ham_adi = CMATRIXList()

    for dof in range(3 * natoms):
        d1ham_adi.append( data_conv.nparray2CMATRIX(per_dof[dof]) )

    return d1ham_adi

def _time_overlap_to_cmatrix(time_overlap):
    return data_conv.nparray2CMATRIX( np.asarray(time_overlap) )

def _nac_vectors_to_dc1_adi(nac_vectors, nstates, natoms ):
    """Pack ES_Result.nac_vectors into Libra's ``dc1_adi``.

    ``nac_vectors`` has shape ``(nstates, nstates, natoms, 3)``; Libra wants one
    CMATRIX(nstates, nstates) per nuclear dof, with ``dof = 3 * atom + xyz``.
    """
    # (nstates, nstates, natoms, 3) -> (3 * natoms, nstates, nstates)
    per_dof = np.asarray(nac_vectors).reshape(nstates, nstates, 3 * natoms)
    per_dof = np.moveaxis(per_dof, -1, 0)

    dc1_adi = CMATRIXList()

    for dof in range(3 * natoms):
        dc1_adi.append( data_conv.nparray2CMATRIX(per_dof[dof]) )

    return dc1_adi

def _time_overlap_to_hvib(  energies, time_overlap, dt):

    nstates = len(energies)
    hvib = _energies_to_ham_adi( energies  )

    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = ( time_overlap[i, j] - time_overlap[j, i] ) / (2.0 * dt)
            hvib.set( i, j, -1j * dij )
            hvib.set( j, i, +1j * dij )

    return hvib

# =============================================================================
# Build Libra callback result
# =============================================================================

def _es_result_to_libra(result: ES_Result, request: ES_Request, natoms: int, dt: float) -> tmp:

    nstates = request.n_total

    obj = tmp()

    # Energies
    obj.ham_adi = _energies_to_ham_adi( result.H_el )

    # Identity adiabatic transformation
    obj.basis_transform = CMATRIX( nstates, nstates )

    for i in range(nstates):
        obj.basis_transform.set( i, i, 1.0 + 0.0j )

    # Gradients
    if result.gradients is not None:
        obj.d1ham_adi = _gradients_to_d1ham_adi( result.gradients, nstates, natoms )

    # Derivative couplings
    if result.nac_vectors is not None:
        obj.dc1_adi = _nac_vectors_to_dc1_adi( result.nac_vectors, nstates, natoms )

    # Time overlaps
    if result.time_overlap is not None:
        obj.time_overlap_adi = _time_overlap_to_cmatrix( result.time_overlap )
        obj.hvib_adi = _time_overlap_to_hvib( result.H_el, result.time_overlap, dt )

    elif request.time_overlap:
        # Overlaps were asked for, but this is the first geometry of the
        # trajectory: there is no previous state to overlap against
        # (ES_Strategy.compute_result only computes one when
        # get_previous_state() is not None).
        #
        # S(t0, t0) is the identity -- no time has elapsed, so every state
        # overlaps perfectly with itself and no coupling has accumulated yet.
        # A zero matrix would instead assert that consecutive electronic states
        # are mutually orthogonal, which is both physically false and rejected
        # outright by nac_npi: validate_npi_input requires |S^T S - I| < 1e-6
        # (src/calculators/NPI.cpp:161), so nac_update_method=2 with
        # nac_algo=NPI would raise on the very first step.
        obj.time_overlap_adi = CMATRIX( nstates, nstates )

        for i in range(nstates):
            obj.time_overlap_adi.set( i, i, 1.0 + 0.0j )

        obj.hvib_adi = _energies_to_ham_adi( result.H_el )

    # Otherwise time-overlaps were never requested, so no attribute is set at
    # all -- the same convention used above for d1ham_adi and dc1_adi.
    # nHamiltonian copies only the attributes that are present (the hasattr
    # gate in nHamiltonian_compute_adiabatic.cpp:544-680), so leaving it off
    # means "this callback has nothing to say about time-overlaps" rather than
    # overwriting whatever the Hamiltonian already holds.

    return obj


# =============================================================================
# Main Libra callback
# =============================================================================

def pyscf_compute_adi(q,params,full_id):

    # 1. Libra input -> generic ES input
    itraj = _get_trajectory_index(full_id)
    geometry = _q_to_geometry(q, itraj, params["atom_labels"] )
    request = _params_to_request(params, itraj)

    # 2. Get this trajectory's persistent ES strategy.
    #
    # ES_Strategy tracks its own previous-geometry snapshot internally
    # (get_previous_state/snapshot_state), so the *same* instance must be reused
    # across calls for a given trajectory or time-overlaps and NACs are never
    # available. Instances are therefore memoized per trajectory index -- the
    # same {itraj: ...} convention the mopac, dftbplus, and cp2k interfaces use
    # for their own per-trajectory caches.

    strategies = params.setdefault("es_strategies", {})

    if itraj not in strategies:
        strategy_spec = params.get(  "strategy_factory",params.get("es_strategy") )

        if strategy_spec is None:
            raise KeyError("Missing strategy specification: expected 'strategy_factory' or 'es_strategy'.")

        if callable(strategy_spec):
            strategies[itraj] = strategy_spec()
        elif hasattr(strategy_spec, "copy"):
            strategies[itraj] = strategy_spec.copy()
        elif hasattr(strategy_spec, "clone"):
            strategies[itraj] = strategy_spec.clone()
        else:
            strategies[itraj] = strategy_spec

    current = strategies[itraj]

    # 3. Run ES calculation.
    result = strategies[itraj].compute_result(geometry, request)

    # 4. Generic ES result -> Libra result
    obj = _es_result_to_libra(
        result,
        request,
        natoms=len(params["atom_labels"]),
        dt=float(params.get("dt", 41.0)),
    )

    return obj

if __name__ == "__main__":
    import os
    import sys

    from liblibra_core import Random, Universe
    import libra_py
    from libra_py import LoadPT
    import libra_py.dynamics.tsh.compute as tsh_dynamics
    from libra_py.packages.pyscf.implementations.casscf import CASSCF

    TSH_FSSH = 0                    # -1 adiabatic (no hops), 0 FSSH, 1 GFSH, 2 MSSH

    # momenta_rescaling_algo -- what direction the momenta are rescaled along
    RESCALE_NONE = 0                # don't rescale
    RESCALE_ALONG_NAC = 201         # along derivative coupling vectors, reverse on frustrated hops
    RESCALE_ALONG_GRAD_DIFF = 211   # along the difference of state-specific forces

    # nac_update_method -- where the NACs come from
    NAC_FROM_DC1 = 1                # contract dc1_adi with p/M  (needs nacv from the ES code)
    NAC_FROM_TIME_OVERLAP = 2       # finite-difference the time-overlap matrix

    # nac_algo -- only in effect when nac_update_method == NAC_FROM_TIME_OVERLAP
    NAC_ALGO_EXTERNAL = -1          # NACs come from somewhere else
    NAC_ALGO_NPI = 1                # Meek & Levine norm-preserving interpolation (0 = HST)

    # force_method
    FORCE_NONE = 0                  # no forces at all (NBRA)
    FORCE_STATE_SPECIFIC = 1        # the active state's force (TSH); 2 = Ehrenfest

    # time_overlap_method -- 0 means "the model supplies time_overlap_adi", which
    # is exactly what this adapter does in _es_result_to_libra.
    TIME_OVERLAP_FROM_MODEL = 0

    # state_tracking_algo -- how the adiabatic states are re-projected step to
    # step. -1 (local diabatization) is the Libra default, and it inverts the
    # time-overlap matrix on every trajectory every step (src/dyn/dyn_ham.cpp:392):
    # a config that never asks the ES code for time-overlaps hands it a zero
    # matrix here, and the run exits with "Problem inverting time-overlap matrix".
    LD_STATE_TRACKING = -1

    REP_ADIABATIC = 1               # rep_* and elec_params["rep"]; 0 = diabatic
    ELEC_INIT_ON_ISTATE = 0         # all trajectories on istate, identical phases
    NUCL_INIT_EXACT = 0             # coords and momenta set exactly to the given values

    # -------------------------------------------------------------------------
    # Shared construction: the parts every config builds the same way
    # -------------------------------------------------------------------------

    def copy_strategies(strategy, ntraj):
        """The ES template, copied once per trajectory.

        Each trajectory needs its own instance: ES_Strategy carries a
        previous-geometry snapshot internally, so sharing one would destroy the
        time-overlaps. Copy while the template is still fresh -- CASSCF.copy() is
        a deepcopy, cheap before the first calculation and expensive after it,
        once mol/mf/mc are attached.
        """
        return {itraj: strategy.copy() for itraj in range(ntraj)}


    def nucl_params_from_geometry(geom, force_constant=0.01):
        """nucl_params straight out of the geometry -- masses included.

        The masses come from Libra's own periodic table (libra_py/elements.dat),
        which LoadPT.Load_PT already converts from Dalton to atomic units
        (LoadPT.py:25). So the atom labels are the only statement of which atoms
        these are; nothing restates their masses by hand.
        """
        U = Universe()
        LoadPT.Load_PT(U, os.path.join(os.path.dirname(libra_py.__file__), "elements.dat"), 0)

        masses = []
        for symbol in geom.atom_labels:
            masses += [U.Get_Element(symbol).Elt_mass] * 3

        ndof = 3 * len(geom.atom_labels)

        # ---- where a sampler plugs in ---------------------------------------
        # NUCL_INIT_EXACT (0) puts every trajectory at exactly this geometry with
        # zero momenta. Libra samples for you with init_type 1/2/3 (widths built
        # from force_constant and mass) or 4 (explicit q_width/p_width lists) --
        # src/dyn/dyn_variables_nuclear.cpp:265-295.
        #
        # For samples drawn OUTSIDE Libra -- a Wigner distribution, say -- use
        # init_type 5, which copies q_init/p_init straight into the dynamical
        # variables with no sampling of its own (dyn_variables_nuclear.cpp:304-313).
        # Both are ndof lists of ntraj floats: the TRANSPOSE of the (ntraj, ndof)
        # arrays a sampler normally hands back.
        #
        #     q, p, _ = wigner_sample(omega, L, masses, R0, ntraj, temperature)
        #     nucl["init_type"] = 5
        #     nucl["q_init"] = q.T.tolist()   # ndof lists of ntraj floats
        #     nucl["p_init"] = p.T.tolist()
        #
        # "q"/"p" below stay as the means: init_type 5 ignores them, but the
        # length of "q" is what fixes ndof for the checks upstream.
        #
        # libra_py.wigner already builds these: generate_wigner_ics() returns one
        # {"traj", "q", "p"} dict per trajectory, and prepare_wigner_from_modes /
        # generate_wigner_from_hessian produce the modes it needs -- from the same
        # masses this function looks up.

        return {
            "ndof": ndof,
            "init_type": NUCL_INIT_EXACT,
            "q": np.asarray(geom.coords_bohr, dtype=np.float64).reshape(-1).tolist(),
            "p": [0.0] * ndof,
            "mass": masses,
            "force_constant": [force_constant] * ndof,
        }


    # -------------------------------------------------------------------------
    # Scenario configurations
    #
    # Each returns the five arguments of generic_recipe, in its argument order:
    #
    #     dyn_params, compute_model, model_params, elec_params, nucl_params
    # -------------------------------------------------------------------------

    def config_fssh_nacv(strategy, geom, ntraj=2, istate=1, dt=41.0, nsteps=2):
        """FSSH with derivative couplings: rescale momenta along the NAC vector.

        Consistent by construction: rescaling along d_ij needs dc1_adi, which the
        ES code only returns when nacv=True; and force_method=1 only ever needs
        the active root's gradient.
        """
        # What these algorithms demand of the ES code -- decided here, with them.
        gradient_state = "active"          # for FORCE_STATE_SPECIFIC
        nacv = True                        # for RESCALE_ALONG_NAC and NAC_FROM_DC1
        time_overlap = True                # NOT for the NACs -- those come from
                                           # dc1_adi -- but for state_tracking_algo
                                           # below, which inverts S every step

        strategies = copy_strategies(strategy, ntraj)
        nucl_params = nucl_params_from_geometry(geom)

        dyn_params = {
            # --- surface hopping ---
            "tsh_method": TSH_FSSH,
            "hop_acceptance_algo": 0,           # based on adiabatic energy
            "momenta_rescaling_algo": RESCALE_ALONG_NAC,
            "use_Jasper_Truhlar_criterion": 1,  # only in effect for algo 201

            # --- Hamiltonian / couplings ---
            "isNBRA": 0,
            "is_nbra": 0,                       # generic_recipe reads this spelling too
            "rep_tdse": REP_ADIABATIC,
            "rep_sh": REP_ADIABATIC,
            "rep_force": REP_ADIABATIC,
            "force_method": FORCE_STATE_SPECIFIC,
            "ham_update_method": 2,             # the model returns the adiabatic Ham directly
            "ham_transform_method": 0,          # no further transformation
            "time_overlap_method": TIME_OVERLAP_FROM_MODEL,
            "state_tracking_algo": LD_STATE_TRACKING,
            "nac_update_method": NAC_FROM_DC1,
            "nac_algo": NAC_ALGO_EXTERNAL,      # unused when nac_update_method == 1
            "hvib_update_method": 1,            # Hvib = Ham - i*hbar*NAC

            # --- integration ---
            "dt": dt,
            "nsteps": nsteps,
            "progress_frequency": 1.0,          # print_freq = int(progress_frequency
                                                # * nsteps) (tsh/save.py:1027), which
                                                # is a division by zero for any run
                                                # shorter than 1/progress_frequency
                                                # steps -- the 0.1 default needs 10+
            "ntraj": len(strategies),
            "quantum_dofs": list(range(nucl_params["ndof"])),

            # --- output ---
            "prefix": "run_fssh_nacv",
            "prefix2": "run_fssh_nacv_aux",
        }
        elec_params = {
            "ndia": strategy.nroots, "nadi": strategy.nroots,
            "init_type": ELEC_INIT_ON_ISTATE,
            "istate": istate,
            "rep": REP_ADIABATIC,
        }
        model_params = {
            "model0": 0,                        # generic_recipe demands it; the callback ignores it
            "atom_labels": geom.atom_labels,
            "dt": dt,
            "nstates": strategy.nroots,         # the template decides the state count

            "gradient_state": gradient_state,
            "nacv": nacv,
            "time_overlap": time_overlap,

            # The per-trajectory strategies, and the initial active state seeded
            # as generic_recipe would from elec_params["istate"]; the driver
            # rewrites act_state every step (dynamics/tsh/compute.py:1082).
            "es_strategies": strategies,
            "act_state": {itraj: istate for itraj in range(ntraj)},
        }
        return dyn_params, pyscf_compute_adi, model_params, elec_params, nucl_params


    def config_nbra(strategy, geom, ntraj=2, istate=0, dt=41.0, nsteps=2):
        """NBRA-style: energies and time-overlaps only, no forces at all.

        The trajectories are propagated on frozen nuclei, so nothing asks the ES
        code for a gradient and the couplings come from the time-overlap matrix.
        See the comment on isNBRA below for why the Hamiltonian is still
        per-trajectory rather than shared.
        """
        gradient_state = None              # for FORCE_NONE
        nacv = False                       # no derivative couplings in the NBRA
        time_overlap = True                # for NAC_FROM_TIME_OVERLAP

        strategies = copy_strategies(strategy, ntraj)
        nucl_params = nucl_params_from_geometry(geom)

        dyn_params = {
            # --- surface hopping ---
            "tsh_method": TSH_FSSH,
            "hop_acceptance_algo": 0,
            "momenta_rescaling_algo": RESCALE_NONE,   # no momenta to rescale in the NBRA
            "use_Jasper_Truhlar_criterion": 0,

            # --- Hamiltonian / couplings ---
            # isNBRA = 1 is what would let ONE Hamiltonian serve every
            # trajectory -- the point of the NBRA. It is off here because it
            # crashes in this build: generic_recipe then allocates a single child
            # (dynamics/tsh/compute.py:1214) while dyn_variables keeps looping
            # over ntraj of them (src/dyn/dyn_variables_electronic.cpp:148), so
            # any ntraj > 1 reads past the end of children[] and segfaults before
            # the first ES call. Set it to 1 only together with ntraj = 1.
            #
            # What the NBRA means for the ES code -- no forces, couplings from
            # time-overlaps alone -- is set below and is independent of the flag.
            "isNBRA": 0,
            "is_nbra": 0,
            "rep_tdse": REP_ADIABATIC,
            "rep_sh": REP_ADIABATIC,
            "rep_force": REP_ADIABATIC,
            "force_method": FORCE_NONE,
            "ham_update_method": 2,
            "ham_transform_method": 0,
            "time_overlap_method": TIME_OVERLAP_FROM_MODEL,
            "state_tracking_algo": LD_STATE_TRACKING,
            "nac_update_method": NAC_FROM_TIME_OVERLAP,
            "nac_algo": NAC_ALGO_NPI,
            "hvib_update_method": 1,

            # --- integration ---
            "dt": dt,
            "nsteps": nsteps,
            "progress_frequency": 1.0,          # print_freq = int(progress_frequency
                                                # * nsteps) (tsh/save.py:1027), which
                                                # is a division by zero for any run
                                                # shorter than 1/progress_frequency
                                                # steps -- the 0.1 default needs 10+
            "ntraj": len(strategies),
            "quantum_dofs": list(range(nucl_params["ndof"])),

            # --- output ---
            "prefix": "run_nbra",
            "prefix2": "run_nbra_aux",
        }
        elec_params = {
            "ndia": strategy.nroots, "nadi": strategy.nroots,
            "init_type": ELEC_INIT_ON_ISTATE,
            "istate": istate,
            "rep": REP_ADIABATIC,
        }
        model_params = {
            "model0": 0,
            "atom_labels": geom.atom_labels,
            "dt": dt,
            "nstates": strategy.nroots,

            "gradient_state": gradient_state,
            "nacv": nacv,
            "time_overlap": time_overlap,

            "es_strategies": strategies,
            "act_state": {itraj: istate for itraj in range(ntraj)},
        }
        return dyn_params, pyscf_compute_adi, model_params, elec_params, nucl_params


    # -------------------------------------------------------------------------
    # The run -- one generic_recipe call, whichever config is selected
    # -------------------------------------------------------------------------

    CONFIGS = { "fssh_nacv": config_fssh_nacv, "nbra": config_nbra }

    name = sys.argv[1] if len(sys.argv) > 1 else "fssh_nacv"
    if name not in CONFIGS:
        sys.exit(f"unknown config {name!r}; pick one of {', '.join(CONFIGS)}")

    # The ES template: 2 electrons in 2 CAS orbitals, 3 roots, HeH+ in STO-3G.
    # Copied once per trajectory by the config; nothing else states the state count.
    strategy = CASSCF(norbcas=2, nelecas=2, nroots=3, basis="sto-3g", charge=1, unit="Bohr")

    # The geometry: HeH+ at its equilibrium bond length, He at the origin.
    # The labels alone fix ndof, the masses, and what the ES code is handed.
    geom = MolecularGeometry(
        atom_labels=("He", "H"),
        coords_bohr=np.array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 1.46379]], dtype=np.float64),
    )

    dyn_params, compute_model, model_params, elec_params, nucl_params = \
        CONFIGS[name](strategy, geom, ntraj=2)

    print(f"config {name}: {''.join(geom.atom_labels)} / {type(strategy).__name__}, "
          f"nstates={elec_params['nadi']}, ntraj={dyn_params['ntraj']}, "
          f"istate={elec_params['istate']}, nsteps={dyn_params['nsteps']}, "
          f"dt={dyn_params['dt']} a.u.")
    print(f"  masses (a.u.): {nucl_params['mass'][::3]}")

    rnd = Random()

    res = tsh_dynamics.generic_recipe(dyn_params, compute_model, model_params,
                                      elec_params, nucl_params, rnd)
