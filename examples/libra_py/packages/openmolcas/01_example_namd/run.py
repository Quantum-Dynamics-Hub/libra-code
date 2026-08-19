import os
import sys
import subprocess
import re
import copy
import numpy as np
from liblibra_core import MATRIX, CMATRIX, CMATRIXList, Py2Cpp_int, Cpp2Py, Random

# Fisrt, we add the location of the library to test to the PYTHON path
import libra_py.packages.dftbplus.methods as dftb
from libra_py.packages.dftbplus.methods import write_dftb_gen, \
            read_dftb_orbital_info, read_overlap_matrix, parse_tagged_file, read_nacv
from libra_py.packages.cp2k import methods as cp2k
import libra_py.citools.ci as ci
from libra_py import data_conv
from libra_py import units
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
from libra_py import initial_conditions

import libra_py.packages.openmolcas.methods as openmolcas

#%matplotlib inline

def do_al3():
    labels = ['Al', 'Al', 'Al']
    coords = [   -1.03293,        1.44509,        0.00000,
                 -2.10052,       -0.80222,        0.00000,
                  0.37950,       -0.60313,        0.00000
             ]

    ndof = len(coords)

    q = [x*units.Angst for x in coords ]

    p = [0.0 for _ in range(ndof)]

    mass = []
    mass_dict = {"Al":26.982 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)

    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")

    nexcitations = 7
    nstates = nexcitations + 1
    return labels, q, p, mass, nat, nexcitations, nstates, ndof



def do_adamantane():
    labels = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 
          'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

    coords = [      9.1785800000,        5.9286800000,        6.7842300000,
           8.2887800000,        6.8273300000,        5.8957000000,
           9.1844900000,        7.7128700000,        5.0000100000,
          10.0683600000,        6.8155300000,        4.1042700000,
          10.9568600000,        5.9168300000,        4.9940500000,
           9.1608400000,        4.1444800000,        5.0000600000,
          10.0624600000,        5.0300200000,        5.8898200000,
           7.3884800000,        5.9404800000,        5.0058800000,
           8.2710700000,        5.0418400000,        4.1102000000,
           9.1667500000,        5.9287100000,        3.2157700000,
           8.5418210000,        5.2981400000,        7.4430300000,
           9.8197100000,        6.5592600000,        7.4387600000,
           7.6533400000,        7.4733510000,        6.5396200000,
           8.5520100000,        8.3737900000,        4.3673800000,
           9.8257000000,        8.3652800000,        5.6327100000,
          10.7122700000,        7.4531000000,        3.4603400000,
          11.6198800000,        6.5472000000,        5.6265800000,
          11.6071300000,        5.2777000000,        4.3570700000,
           9.7891300000,        3.4835800000,        4.3632500000,
           8.5238300000,        3.4920700000,        5.6369700000,
          10.7021300000,        4.3839800000,        6.5295200000,
           6.7297000000,        5.3100700000,        5.6427900000,
           6.7339500000,        6.5795900000,        4.3732500000,
           7.6229200000,        4.4042900000,        3.4705200000,
           9.7951300000,        5.2897900000,        2.5569900000,
           8.5340190000,        6.5676900000,        2.5612300000
    ]

    ndof = len(coords)

    q = [x*units.Angst for x in coords ]

    p = [0.0 for _ in range(ndof)]

    mass = []
    mass_dict = {"C":12.0, "H":1.0 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)

    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")

    nexcitations = 5
    nstates = nexcitations + 1
    return labels, q, p, mass, nat, nexcitations, nstates, ndof



def do_cyclopropanone():
    labels = ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'O']
    coords = [ -0.05118556,  0.16553701,  -0.16798526,
    0.35744328,  0.25865008,  1.31797980,
    1.31544856,  -0.14558136,  0.27790752,
    -0.72277666,  -0.65457267,  -0.44737499,
    -0.05105795,  -0.50171510,  1.99358658,
    0.42993677,  1.26309033,  1.75077019,
    -0.24176175,  1.11030703,  -0.69025457,
    2.42779862,  -0.48800465,  -0.00704308
         ]
    ndof = len(coords)

    q = [x*units.Angst for x in coords ]

    p = [0.0 for _ in range(ndof)]

    mass = []
    mass_dict = {"C":12.0, "H":1.0, "O": 16.0 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)

    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")

    nexcitations = 5
    nstates = nexcitations + 1
    
    return labels, q, p, mass, nat, nexcitations, nstates, ndof


def do_pyrazine():
    labels = ['N', 'N', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H']
    
    coords = [
         1.397389772,  -0.0008423226,  -0.0000282719,
        -1.397375322,   0.0008217057,  -0.0000247962,
         0.7002955864, -1.1432855270,  -0.0000115720,
        -0.7016386487, -1.1424470400,  -0.0000347895,
         0.7016531643,  1.1424276620,  -0.0000302877,
        -0.7002809311,  1.1432631340,  -0.0000160692,
         1.2652815770, -2.0944791720,  -0.0000042173,
        -1.2677524290, -2.0929674790,  -0.0000728923,
         1.2677740730,  2.0929461430,  -0.0000726997,
        -1.2652605130,  2.0944600380,  -0.0000044041
    ]
    
    ndof = len(coords)
    
    q = [x*units.Angst for x in coords ]
    
    p = [0.0 for _ in range(ndof)]
    
    mass = []
    mass_dict = {"C":12.0, "H":1.0, "N": 14.0 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)
    
    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")
    
    
    nexcitations = 5
    nstates = nexcitations + 1
    return labels, q, p, mass, nat, nexcitations, nstates, ndof


def do_cyclobutanone():
    labels = ['O', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']

    coords = [
        1.89203466, -0.00001005,  0.10213634,
       -0.37508690, -1.10363651,  0.01410641,
       -1.46025913, -0.00000405, -0.03142451,
       -0.37509099,  1.10362853,  0.01419800,
        0.69593416, -0.00000447,  0.05578079,
       -0.33350719,  1.75085441, -0.87700831,
       -0.40547773,  1.74520723,  0.90990079,
       -2.06405275,  0.00002296, -0.95243499,
       -2.13881477, -0.00003568,  0.83585874,
       -0.33348542, -1.75074241, -0.87718510,
       -0.40547813, -1.74534221,  0.90971976
    ]

    ndof = len(coords)

    q = [x*units.Angst for x in coords ]

    p = [0.0 for _ in range(ndof)]

    mass = []
    mass_dict = {"C":12.0, "H":1.0, "O": 16.0 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)

    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")
    nexcitations = 4
    nstates = nexcitations + 1
    return labels, q, p, mass, nat, nexcitations, nstates, ndof

def do_furan():
    labels = ['C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'O']
    
    coords = [
        -4.17442176, -0.42652134, -1.10769802,
        -3.89378544, -0.03496919, -2.08932411,
        -4.94077927, -1.49806576, -0.71470450,
        -5.44969324, -2.21007698, -1.36814310,
        -4.17447024, -0.42647113,  1.10773382,
        -3.89386725, -0.03488139,  2.08935445,
        -4.94081052, -1.49803339,  0.71475532,
        -5.44976392, -2.21000728,  1.36820389,
        -3.70250447,  0.23335116,  0.00001327
    ]
    
    ndof = len(coords)
    
    q = [x*units.Angst for x in coords ]
    
    p = [0.0 for _ in range(ndof)]
    
    mass = []
    mass_dict = {"C":12.0, "H":1.0, "O": 16.0 }
    for elt in labels:
        m = mass_dict[elt] * units.amu
        mass.append(m)
        mass.append(m)
        mass.append(m)
    
    print(F"ndof = {ndof}")
    nat = len(labels)
    print(F"nat = {nat}")
    nexcitations = 4
    nstates = nexcitations + 1
    return labels, q, p, mass, nat, nexcitations, nstates, ndof



#labels, q, p, mass, nat, nexcitations, nstates, ndof = do_adamantane()
#labels, q, p, mass, nat, nexcitations, nstates, ndof = do_cyclopropanone()
#labels, q, p, mass, nat, nexcitations, nstates, ndof = do_pyrazine()
#labels, q, p, mass, nat, nexcitations, nstates, ndof = do_cyclobutanone()
#labels, q, p, mass, nat, nexcitations, nstates, ndof = do_furan()
labels, q, p, mass, nat, nexcitations, nstates, ndof = do_al3()

ntraj = 3
DT = 0.5 # fs
ISTATE = 2


############################################################
### 1. Choose the initial conditions: Nuclear and Electronic
############################################################

rng = np.random.default_rng()

q_fluct = rng.random( (ndof, ntraj) )
p_fluct = rng.random( (ndof, ntraj) )

q_init = [ [ q[idof] +  0.1*q_fluct[idof, itraj] for itraj in range(ntraj)] for idof in range(ndof) ]
p_init = [ [ p[idof] + 0.05*p_fluct[idof, itraj] for itraj in range(ntraj)] for idof in range(ndof) ]

# Remove total linear and angular momenta:
q_init_np = np.array(q_init)
p_init_np = np.array(p_init)
mass_np = np.array(mass)

p_init_np = initial_conditions.cleanup_momenta(q_init_np, p_init_np, mass_np)
p_init = p_init_np.tolist()


print(q_init)
print(p_init)

#sys.exit(0)

# Nuclear initial conditions 
nucl_params = { "ndof":ndof,
                "q":q, "p":p,
                "mass":mass,
                "force_constant":[ 0 for _ in range(ndof) ],
                "q_width":[ 0.1 for _ in range(ndof)  ],
                "p_width":[ 0.1 for _ in range(ndof)  ],
                "q_init" : q_init,
                "p_init" : p_init,
                "init_type": 5
              }

# Electronic initial conditions - start on the second adiabatic state (starting from 0 = ground)
##########
istate = ISTATE
##########
istates = [0.0 for i in range(nstates)]
istates[istate] = 1.0     
elec_params = {"verbosity":2, "init_dm_type":0,"ndia":nstates, "nadi":nstates, 
               "rep":1, "init_type":1, "istates":istates, "istate":istate     }


###########################################################
### 2. Dynamics parameters
###########################################################

dyn_general = { "nsteps":5, "ntraj":ntraj, "nstates":nstates,
                "dt":DT*units.fs2au, "num_electronic_substeps":1, "isNBRA":0, "is_nbra":0,
                "progress_frequency":0.5, "which_adi_states":range(nstates), "which_dia_states":range(nstates),
                "mem_output_level":3,
                "properties_to_save":[ "timestep", "time", "q", "p", "f", "Cadi", "Cdia", "Epot_ave", "Ekin_ave", "Etot_ave",
                "states", "se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia"],
                "prefix":"adiabatic_md", "prefix2":"adiabatic_md"
              }

from recipes import fssh2
fssh2.load(dyn_general) # FSSH2


###########################################################
### 3. Model parameters
###########################################################
#os.chdir("/home/alexvakimov/Projects/Project_DFTB")
wd = "al3"

# Create working directory, if doesn't exist
if not os.path.exists(wd):
    os.mkdir(wd)


couplings = "{0 " + F"{nexcitations}" + "}" 

dftb_run_params = {
    "gen_file" : "x1.gen",
    "sk_prefix" : "../mio/FinalSK/",
    "Driver" : "{}",
    "MaxAngularMomentum" : """{ C = "p"
                                H = "s"
                                O = "p"
                     }
                     """,
    "Symmetry" : "Singlet",
    "NrOfExcitations" : nexcitations,
    "StateOfInterest" : 1,
    "WriteSPTransitions" : "Yes",
    "WriteXplusY" : "Yes",
    "WriteXplusYAscii" : "Yes",
    "StateCouplings" : couplings,
    
    "WriteAutotestTag" : "Yes",
    "WriteHS" : "No",
    "WriteEigenvectors" : "Yes",
    "EigenvectorsAsText" : "Yes",
    "PrintForces" : "Yes",
    "Filling" : """Fermi { Temperature [K] = 0.01 }
                """
}


molcas_run_params = {
        "basis"      : "cc-pVDZ",
        "charge"     : 0,
        "spin"       : 2,
        "scf_method" : "uhf",
        "uhf_orbital_set" : "beta",
        "title"      : "al3",

        "nactel"     : "3 0 0",
        "inactive"   : 18,
        "ras2"       : 6,
        "ciroot"     : f"{nstates} {nstates} 1",

        "nac_states" : None,
        "group"      : "NoSym",
        "prweight"   : 1e-6,
        "thre"       : 1e-6,
        "lshift"     : 0.1,            # Level shift for convergence
        "maxiter"    : 500,
    }


#print(dftb_params)
print(molcas_run_params)
#sys.exit(0)


model_params_dftb = {"atom_labels":labels, "timestep":0, 
               "dftb_exe":"/home/alexvakimov/SOFTWARE/dftbplus/_install/bin/dftb+", 
               "dftb_run_params": dftb_run_params,
               "working_directory_prefix":"non_nbra_wd",
               
               "odin_exe": "/home/alexvakimov/SOFTWARE/odin/odin",
               "odin_max_ang_mom" : { "C":2, "H":1, "O":2, "N":2 },
               "orbital_space" : None,
               
               "dt":DT*units.fs2au,
               "nelec_act_space":None,
               "ci_threshold":0.01,
               
               "act_state": { i:istate for i in range(ntraj) },
               "is_first_time" : { i:True for i in range(ntraj) },
               
               "read_forces" : True,
               "read_nacvs": True,
               
               "nstates": nstates,
                "model": 0,
                "model0": 0
              }

model_params_molcas = {
    "atom_labels": labels,
    "timestep": 0,
    "dt": DT * units.fs2au,
    "exe": "/home/alexvakimov/SOFTWARE/OpenMolcas/_build/pymolcas",
    "molcas_run_params": molcas_run_params,
    "working_directory_prefix": "wd_molcas",
    "molcas_input_prefix": "input_",
    "molcas_output_prefix": "output_",
    "verbose": True,
    "do_Lowdin": True,

    # Total Libra dimension after M_S expansion
    "nstates": 8,

    # |C|^2 threshold, shared by all manifolds
    "ci_threshold": 1e-6,

    "act_state": {i: istate for i in range(ntraj)},
    "is_first_time": {i: True for i in range(ntraj)},

    "spin_manifolds": [
        {
            "spin": 2,
            "nroots": 2,
            "ciroot": "2 2 1",
            "gradient_states": "all",
            "nac_pairs": "all",
            "nac_nocsf": False,
        },
        {
            "spin": 4,
            "nroots": 1,
            "ciroot": "1 1 1",
            "gradient_states": "all",
            "nac_pairs": "all",
            "nac_nocsf": False,
        },
    ],

    "model": 0,
    "model0": 0,
}

# The ordering is:
"""
obj.spin_labels == [
    (2, 0,  1),   # doublet root 0, Ms = +1/2
    (2, 0, -1),   # doublet root 0, Ms = -1/2
    (2, 1,  1),   # doublet root 1, Ms = +1/2
    (2, 1, -1),   # doublet root 1, Ms = -1/2

    (4, 0,  3),   # quartet root 0, Ms = +3/2
    (4, 0,  1),   # quartet root 0, Ms = +1/2
    (4, 0, -1),   # quartet root 0, Ms = -1/2
    (4, 0, -3),   # quartet root 0, Ms = -3/2
]
"""

#model_params = model_params_dftb
model_params = model_params_molcas


dyn_params = dict(dyn_general)
pref = "FSSH2_"
dyn_params.update({ "prefix":pref, "prefix2":pref })
print(F"Computing {pref}")    


rnd = Random()
#res = tsh_dynamics.generic_recipe(dyn_params, dftb.dftb_compute_adi, model_params, elec_params, nucl_params, rnd)
res = tsh_dynamics.generic_recipe(dyn_params, openmolcas.molcas_compute_adi, model_params, elec_params, nucl_params, rnd)


############################################################
### 5. Plot the results - preliminary plotting
############################################################
pref = "FSSH2_"
NSTATES = nstates

NTRAJ = dyn_general["ntraj"]
plot_params = { "prefix":pref, "filename":"mem_data.hdf", "output_level":3,
                "which_trajectories":list(range(NTRAJ)), "which_dofs":[0], "which_adi_states":list(range(NSTATES)), 
                "which_dia_states":list(range(NSTATES)), 
                "frameon":True, "linewidth":3, "dpi":300,
                "axes_label_fontsize":(8,8), "legend_fontsize":8, "axes_fontsize":(8,8), "title_fontsize":8,
                "what_to_plot":["coordinates", "momenta",  "forces", "energies", "phase_space", "se_pop_adi",
                                "se_pop_dia", "sh_pop_adi", "sh_pop_dia" ], 
                "which_energies":["potential", "kinetic", "total"],
                "save_figures":1, "do_show":0, "no_label":1
              }
tsh_dynamics_plot.plot_dynamics(plot_params)


import h5py
import matplotlib.pyplot as plt
from libra_py import units

#%matplotlib inline

with h5py.File("FSSH2_/mem_data.hdf", 'r') as f:
    t = f["time/data"][:]
    energy = f["Etot_ave/data"][:] - f["Etot_ave/data"][0]
    epot = f["Epot_ave/data"][:] - f["Epot_ave/data"][0]
    ekin = f["Ekin_ave/data"][:] - f["Ekin_ave/data"][0]

    plt.plot(t * units.au2fs, energy * units.au2ev, label="Total")
    plt.plot(t * units.au2fs, epot * units.au2ev, label="Pot")
    plt.plot(t * units.au2fs, ekin * units.au2ev, label="Kin")

    plt.legend()
    #plt.show()
    plt.savefig("FSSH2_/dE.png")
