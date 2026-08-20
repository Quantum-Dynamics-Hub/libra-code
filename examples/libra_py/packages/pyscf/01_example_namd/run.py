"""Mixed-multiplicity CASSCF/FSSH dynamics for an Al3 cluster."""

import os

import numpy as np

from liblibra_core import Random
from libra_py import initial_conditions, units
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.methods import pyscf_compute_adi


def do_al3():
    """Return the Al3 geometry and nuclear initial-condition parameters."""
    labels = ["Al", "Al", "Al"]
    coords_angstrom = [
        -1.03293,  1.44509, 0.00000,
        -2.10052, -0.80222, 0.00000,
         0.37950, -0.60313, 0.00000,
    ]
    q = [value * units.Angst for value in coords_angstrom]
    p = [0.0] * len(q)
    mass = []
    for _ in labels:
        mass.extend([26.982 * units.amu] * 3)
    return labels, q, p, mass


labels, q, p, mass = do_al3()
ndof = len(q)
nat = len(labels)
ntraj = int(os.environ.get("LIBRA_NTRAJ", "3"))
nstates = 8
DT = 0.5  # fs
ISTATE = 2
nsteps = int(os.environ.get("LIBRA_NSTEPS", "5"))

print(f"ndof = {ndof}")
print(f"nat = {nat}")


############################################################
### 1. Choose the initial conditions: Nuclear and Electronic
############################################################

rng = np.random.default_rng()
q_fluct = rng.random((ndof, ntraj))
p_fluct = rng.random((ndof, ntraj))

q_init = [
    [q[idof] + 0.1 * q_fluct[idof, itraj] for itraj in range(ntraj)]
    for idof in range(ndof)
]
p_init = [
    [p[idof] + 0.05 * p_fluct[idof, itraj] for itraj in range(ntraj)]
    for idof in range(ndof)
]

p_init = initial_conditions.cleanup_momenta(
    np.asarray(q_init), np.asarray(p_init), np.asarray(mass)
).tolist()

nucl_params = {
    "ndof": ndof,
    "q": q,
    "p": p,
    "mass": mass,
    "force_constant": [0.0] * ndof,
    "q_width": [0.1] * ndof,
    "p_width": [0.1] * ndof,
    "q_init": q_init,
    "p_init": p_init,
    "init_type": 5,
}

istates = [0.0] * nstates
istates[ISTATE] = 1.0
elec_params = {
    "verbosity": 2,
    "init_dm_type": 0,
    "ndia": nstates,
    "nadi": nstates,
    "rep": 1,
    "init_type": 1,
    "istates": istates,
    "istate": ISTATE,
}


###########################################################
### 2. Dynamics parameters
###########################################################

dyn_general = {
    "nsteps": nsteps,
    "ntraj": ntraj,
    "nstates": nstates,
    "dt": DT * units.fs2au,
    "num_electronic_substeps": 1,
    "isNBRA": 0,
    "is_nbra": 0,
    "progress_frequency": max(0.5, 1.0 / nsteps),
    "which_adi_states": range(nstates),
    "which_dia_states": range(nstates),
    "mem_output_level": 3,
    "properties_to_save": [
        "timestep", "time", "q", "p", "f", "Cadi", "Cdia",
        "Epot_ave", "Ekin_ave", "Etot_ave", "states",
        "se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia",
    ],
    "prefix": "adiabatic_md",
    "prefix2": "adiabatic_md",
}

from recipes import fssh2

fssh2.load(dyn_general)
# The recipe normally rescales along explicit NAC vectors. PySCF's
# spin-constrained state-averaged CASSCF NAC response is not available, while
# all state gradients are, so use the gradient-difference direction instead.
dyn_general["momenta_rescaling_algo"] = 211


###########################################################
### 3. PySCF model parameters
###########################################################

doublet = CASSCF(
    norbcas=3,
    nelecas=(2, 1),
    nroots=2,
    basis="sto-3g",
    charge=0,
    unit="Bohr",
    spin_multiplicity=2,
)
quartet = CASSCF(
    norbcas=3,
    nelecas=(3, 0),
    nroots=1,
    basis="sto-3g",
    charge=0,
    unit="Bohr",
    spin_multiplicity=4,
)

model_params = {
    "atom_labels": labels,
    "dt": DT * units.fs2au,
    "nstates": nstates,
    "energy_zero": -716.7,
    "gradient_state": "all",
    "nacv": False,
    "time_overlap": True,
    "spin_manifolds": [
        {
            "spin": 2,
            "nroots": 2,
            "es_strategy": doublet,
        },
        {
            "spin": 4,
            "nroots": 1,
            "es_strategy": quartet,
        },
    ],
    "model": 0,
    "model0": 0,
}

# State ordering:
#   (2, root 0, Ms=+1/2), (2, root 0, Ms=-1/2),
#   (2, root 1, Ms=+1/2), (2, root 1, Ms=-1/2),
#   (4, root 0, Ms=+3/2), (4, root 0, Ms=+1/2),
#   (4, root 0, Ms=-1/2), (4, root 0, Ms=-3/2).
# Without spin-orbit coupling, different multiplicities and different Ms
# components have zero Hamiltonian, overlap, and derivative-coupling blocks.


###########################################################
### 4. Run dynamics
###########################################################

pref = "FSSH2_"
dyn_params = dict(dyn_general)
dyn_params.update({"prefix": pref, "prefix2": pref})
print(f"Computing {pref}")

res = tsh_dynamics.generic_recipe(
    dyn_params,
    pyscf_compute_adi,
    model_params,
    elec_params,
    nucl_params,
    Random(),
)


############################################################
### 5. Plot the results
############################################################

plot_params = {
    "prefix": pref,
    "filename": "mem_data.hdf",
    "output_level": 3,
    "which_trajectories": list(range(ntraj)),
    "which_dofs": [0],
    "which_adi_states": list(range(nstates)),
    "which_dia_states": list(range(nstates)),
    "frameon": True,
    "linewidth": 3,
    "dpi": 300,
    "axes_label_fontsize": (8, 8),
    "legend_fontsize": 8,
    "axes_fontsize": (8, 8),
    "title_fontsize": 8,
    "what_to_plot": [
        "coordinates", "momenta", "forces", "energies", "phase_space",
        "se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia",
    ],
    "which_energies": ["potential", "kinetic", "total"],
    "save_figures": 1,
    "do_show": 0,
    "no_label": 1,
}
tsh_dynamics_plot.plot_dynamics(plot_params)
