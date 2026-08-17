import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import libra_py.packages.mopac.methods as mopac
import libra_py.packages.cp2k.methods as cp2k
import libra_py.units as units
import libra_py.citools.ci as ci

from liblibra_core import Py2Cpp_int

example_dir = Path(__file__).resolve().parent
output_dir = example_dir / "output"
output_dir.mkdir(exist_ok=True)
wd = str(output_dir / "workflow_wd_standard")

trajectory = example_dir / "TiO2-aligned.xyz"


def read_frame(step):
    """Read a frame while containing files from older Libra installations."""
    try:
        return cp2k.read_trajectory_xyz_file(trajectory, step, output_dir)
    except TypeError:
        previous_directory = Path.cwd()
        try:
            os.chdir(output_dir)
            return cp2k.read_trajectory_xyz_file(str(trajectory), step)
        finally:
            os.chdir(previous_directory)


labels, q = read_frame(0)

params = {"atom_labels":labels, 
          "timestep": 0,
          "exe":"/home/alexvakimov/SOFTWARE/mopac/_build/mopac", 
          "mopac_run_params":"INDO C.I.=(6,3) CHARGE=-1 RELSCF=0.000001 ALLVEC WRTCONF=0.10 WRTCI=4",
          "multiplicity":2,
          "spin_projection":0.5,
          "nelec_act_space":7,
          "working_directory_prefix":wd,
          "mopac_input_prefix":"tio2_", 
          "mopac_output_prefix":"output_", 
          "nstates":4,
          "dt":1.0*units.fs2au,
          "do_Lowdin":True,
               
          "is_first_time":{0:True},
          "act_state":{0:1},
         }

# For 1 trajectory
print(params)

# Emulates 1 trajectory
full_id = Py2Cpp_int([0, 0])

res = F"{wd}_itraj0"

# Do the first 5 steps 
for i in range(5):
    labels, q = read_frame(i)
    params["timestep"] = i
    
    obj = mopac.mopac_compute_adi(q, params, full_id)        
    obj.ham_adi.show_matrix(F"{res}/ham_adi_{i}.txt")
    obj.hvib_adi.show_matrix(F"{res}/hvib_adi_{i}.txt")
    obj.time_overlap_adi.real().show_matrix(F"{res}/st_adi_{i}.txt")
    obj.overlap_adi.real().show_matrix(F"{res}/s_adi_{i}.txt")
