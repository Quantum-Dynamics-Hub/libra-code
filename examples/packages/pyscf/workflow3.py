from liblibra_core import Random
from libra_py.dynamics.tsh.recipes import fssh2_v_plus as fssh2
import libra_py.dynamics.tsh.compute as tsh_dynamics

# Import the Libra-provided universal adapter
from libra_py.packages.compute_adi_adapter import es_compute_adi

# Import the user's specific ES implementation
from my_quantum_code import CASSCF
from src.libra_py import units

# 0. Define the molecule equilibrum geometry and trajectory count
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
mass = [mass_dict[label] * units.amu for label in labels]

# 1. Setup the ES engine (a list of instances for each trajectory)
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

# 2. Setup the Dynamics Params
dyn_params = {
    "recipe": fssh2,
    "dt": 0.1,
    "ntraj": ntraj,
    "nsteps": 100,
    "nstates": nstates,
    "which_adi_states": list(range(nstates)),
    "which_dia_states": list(range(nstates)),
    "num_electronic_substeps": 1,
    "mem_output_level": 3,
    "prefix": "test",
}

# 3. Setup Params
nucl_params = {
    "ndof": ndof,
    "q": q,
    "p": p,
    "mass": mass,
    "init_type": 3,
}
elec_params = {...}

model_params = {
    "model0": 0, 
    "atom_labels": ["Li", "H"],
    "nstates": 4,
    "nat": 2,
    "es_strategy": es_engines, # Pass the ABC instances
}

# 4. Run Dynamics
res = tsh_dynamics.generic_recipe(
    dyn_params,
    es_compute_adi,  # Pass the Libra internal adapter
    model_params,
    elec_params,
    nucl_params,
    Random(),
)