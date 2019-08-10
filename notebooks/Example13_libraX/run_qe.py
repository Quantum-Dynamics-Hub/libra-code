import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import workflows.nbra



user = 1  # 0 - Alexey, 1 - Ekadashi

################ System-specific settings ########################
if user==0:
    # For Alexey
    libra_bin_path = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code/_build/src" # set the path name to the source files in libracode
    libra_x_path = "/user/alexeyak/Programming/Libra-X/src"

elif user==1:
    # For Ekadashi
    libra_bin_path = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build/src"
    libra_x_path = "/projects/academic/alexeyak/ekadashi/Libra-X/src"


os.environ["src_path"] = libra_x_path   # Path to the source code
sys.path.insert(1,os.environ["src_path"])    # Path to the source code


import main # import main module of the libra-QE-interface code
import defaults


########## Setup all manual parameters here ####################

params = {}

defaults.set_defaults(params, "QE")

params["nproc"] = 12              # the number of processors
params["dt_nucl"]=20.0  # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["Nsnaps"]=5    # the number of MD rounds
params["Nsteps"]=2      # the number of MD steps per snap
params["Ncool"] = -1
params["nspin"] = 2
params["electronic_smearing"] = 0.01 # Electronic smearing used in Fermi population calculation
params["nconfig"] = 1
params["el_mts"] = 1
params["num_SH_traj"] = 1
params["scf_itr"] = 10  # Number of SCF steps in each fractional occupation update
params["max_iteration"] = 30 # Maximum number of fractional Fermi cycle
params["non-orth"] = 1  # = 1 when MOs are non-orthogonal, = 0 when calculated in orthogonal MO basis
params["print_S_mat"] = 0 # 1 if S-matrix printing required, 0 if not required
params["smat_inc"] = 0 # 1 Including overlap matrix (S), 0 when overlap matrix (S) not included in el propagation

params["MD_type"] = 1  # 1 NVT ensamble, 0 NVE ensamble
# Thermostat parameters
params["therm"] = Thermostat({"thermostat_type":"Nose-Hoover","nu_therm":0.001,"Temperature":300.0,"NHC_size":3})
params["Temperature"] = params["therm"].Temperature
#params["Temperature"] = 300.0
#params["nu_therm"] = 0.001
#params["NHC_size"] = 3
#params["thermostat_type"] = "Nose-Hoover"
params["sigma_pos"] = 0.01  #Displace atomic position randomly
params["Nstart"] = 0

########### Now start actual calculations ###########################

#params["num_MO"] = 3  # number of MO basis used in constructing electronic wavefunction
params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1) ] 
params["excitations_init"] = [1]
params["HOMO"] = 0
params["min_shift"] = 0
params["max_shift"] = 1 
for i in range(0,len(params["excitations"])):
    params["qe_inp0%i" %i] = "x%i.scf.in" %i    # initial input file
    params["qe_inp%i" %i] = "x%i.scf_wrk.in" %i # working input file 
    params["qe_out%i" %i] = "x%i.scf.out" %i    # output file


main.main(params)  # run actual calculations


