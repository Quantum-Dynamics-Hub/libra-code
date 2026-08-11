from liblibra_core import Random
from libra_py import units
import libra_py.dynamics.tsh.compute as tsh_dynamics
from libra_py.dynamics.tsh.recipes import fssh2_v_plus as fssh2

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.methods import es_compute_adi



# Molecule and trajectory count


labels = ["Li", "H"]

coords_bohr = [
    0.000, 0.000, 0.0000,
    0.000, 0.000, 1.0000,
]

ndof = len(coords_bohr) * 3
nat = len(labels)

nstates = 4
ntraj = 2

q = [x * units.Bohr for x in coords_bohr]
p = [0.0 for _ in range(ndof)]

mass_dict = {"Li": 7.0, "H": 1.0}
mass = []
for label in labels:
    m = mass_dict[label] * units.amu
    mass.extend([m, m, m])

# Nuclear and electronic initial conditions

nucl_params = {
    "ndof": ndof,
    "q": q,
    "p": p,
    "mass": mass,
    "force_constant": [0.0 for _ in range(ndof)],
    "q_width": [0.25 for _ in range(ndof)],
    "p_width": [0.1 for _ in range(ndof)],
    "init_type": 4,
}

istate = 2
istates = [0.0 for _ in range(nstates)]
istates[istate] = 1.0

elec_params = {
#...
}

# Electronic-structure engine instantiation
# Each trajectory needs its own stateful engine for cached ES history.

es_engines = [
    CASSCF(
        norbcas=4,
        nelecas=2,
        nroots=nstates,
        basis="sto-3g",
        charge=0,
    )
    for _ in range(ntraj)
]

# Dynamics parameters

dyn_params = {
    "nsteps": 5,
    "ntraj": ntraj,
    "nstates": nstates,
#...
}

fssh2.load(dyn_params)

# Model params passed to compute_adi

model_params = {
    "model": 0,
    "model0": 0,
    "atom_labels": labels,
    "nstates": nstates,
    "dt": dyn_params["dt"],
    "es_strategy": es_engines,
}

# Optional equilibration stage 

# Run NAMD

rnd = Random()

res = tsh_dynamics.generic_recipe(
    dyn_params,
    es_compute_adi,
    model_params,
    elec_params,
    nucl_params,
    rnd,
)
