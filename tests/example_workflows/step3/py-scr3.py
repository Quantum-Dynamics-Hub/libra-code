from PYXAID import *
import os

#############################################################################################
# Input section: Here everything can be defined in programable way, not just in strict format
#############################################################################################

params = {}

# Define general control parameters (file names, directories, etc.)
# Path to Hamiltonians
# These paths must direct to the folder that contains the results of
# the step2 calculations (Ham_ and (optinally) Hprime_ files) and give
# the prefixes and suffixes of the files to read in
rt = "/projects/academic/alexeyak/brendan/official_workflows_libra/libra-code/tests/example_workflows/step3" 
params["Ham_re_prefix"] = rt+"/res/0_Ham_"
params["Ham_re_suffix"] = "_re"
params["Ham_im_prefix"] = rt+"/res/0_Ham_"
params["Ham_im_suffix"] = "_im"
#params["Hprime_x_prefix"] = rt + "/res/0_Hprime_"
#params["Hprime_x_suffix"] = "x_re"
#params["Hprime_y_prefix"] = rt + "/res/0_Hprime_"
#params["Hprime_y_suffix"] = "y_re"
#params["Hprime_z_prefix"] = rt + "/res/0_Hprime_"
#params["Hprime_z_suffix"] = "z_re"
params["energy_units"] = "Ry"                # This specifies the units of the Hamiltonian matrix elements as they
                                             # are written in Ham_ files. Possible values: "Ry", "eV"

# Set up other simulation parameters:
# Files and directories (apart from the Ham_ and Hprime_)
params["scratch_dir"] =  os.getcwd()+"/out"  # Hey! : you need to create this folder in the current directory
                                             # This is were all (may be too many) output files will be written
params["read_couplings"] = "batch"           # How to read all input (Ham_ and Hprime_) files. Possible values:
                                             # "batch", "online"

# Simulation type
params["runtype"] = "namd"                   # Type of calculation to perform. Possible values:
                                             # "namd" - to do NA-MD calculations, "no-namd"(or any other) - to
                                             # perform only pre-processing steps - this will create the files with
                                             # the energies of basis states and will output some useful information,
                                             # it may be particularly helpful for preparing your input
params["decoherence"] = 1                    # Do you want to include decoherence via DISH? Possible values:
                                             # 0 - no, 1 - yes
params["is_field"] = 0                       # Do you want to include laser excitation via explicit light-matter
                                             # interaction Hamiltonian? Possible values: 0 - no, 1 - yes

# Integrator parameters
params["elec_dt"] = 1.0                      # Electronic integration time step, fs
params["nucl_dt"] = 1.0                      # Nuclear integration time step, fs (this parameter comes from 
                                             # you x.md.in file)
params["integrator"] = 0                     # Integrator to solve TD-SE. Possible values: 0, 10,11, 2

# NA-MD trajectory and SH control 
params["namdtime"] = 75                      # Trajectory time, fs
params["num_sh_traj"] = 1000                  # Number of stochastic realizations for each initial condition
params["boltz_flag"] = 1                     # Boltzmann flag (set to 1 anyways)
params["Temp"] = 300.0                       # Temperature of the system
params["alp_bet"] = 0                        # How to treat spin. Possible values: 0 - alpha and beta spins are not
                                             # coupled to each other, 1 - don't care about spins, only orbitals matter

params["debug_flag"] = 0                     # If you want extra output. Possible values: 0, 1, 2, ...
                                             # as the number increases the amount of the output increases too
                                             # Be carefull - it may result in a huge output!

# Parameters of the field (if it is included)
params["field_dir"] = "xyz"                 # Direction of the field. Possible values: "x","y","z","xy","xz","yz","xyz"
params["field_protocol"] = 1                # Envelope function. Possible values: 1 - step function, 2 - saw-tooth
params["field_Tm"] = 25.0                   # Middle of the time interval during which the field is active
params["field_T"] = 25.0                    # The period (duration) of the field pulse
params["field_freq"] = 3.0                  # The frequency of the field radiation = energy of the photons
params["field_freq_units"] = "eV"           # Units of the above quantity. Possible values: "eV", "nm","1/fs","rad/fs"
params["field_fluence"] = 1.0               # Defines the light radiation intensity (fluence), mJ/cm^2



# Define states:
# Example of indexing convention with Nmin = 5, HOMO = 5, Nmax = 8
# the orbitals indices are consistent with QE (e.g. PP or DOS) indexing, which starts from 1
# [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] - all computed orbitals
# [ 1, 2, 3, 4, 5, 6]                     - occupied orbitals
#                   [ 7, 8, 9, 10, 11] - unoccupied orbitals
#              [5, 6, 7, 8]            - active space


# Set active space and the basis states
params["active_space"] = range(1,22)

params["states"] = []
params["states"].append(["GS",[10,-10,11,-11],0.00])  
params["states"].append(["S1",[10,-10,11,-12],0.00])  
params["states"].append(["S1",[10,-10,11,-13],0.00])  
params["states"].append(["S1",[10,-10,11,-14],0.00])
params["states"].append(["S1",[10,-10,11,-15],0.00])

# Initial conditions
nmicrost = len(params["states"])
ic = []
i = 0
while i<10:
    j = 0
    while j<nmicrost:
        ic.append([2*i,j])
        j = j + 1
    i = i + 1

params["iconds"] = ic

# Scale NACs
# Below is a small snippet that scales all NACs by 1000.0 (so no practical effect)
nmicrost = len(params["states"]) # number of (micro)states
params["nac_scale"] = []
for i in range(0,nmicrost):
    for j in range(0,nmicrost):
        params["nac_scale"].append([i,j,1000.0]) 



#############################################################################################
# Execution section: Here we actually start the NA-MD calculations and the analysis
#############################################################################################

############ Run calculations ######################
print params                   # print out all simulation parameters first
pyxaid_core.info().version()
pyxaid_core.namd(params)


########### Below we will be using the average.py module ########
# Note: If you want to re-run averaging calculations - just comment out the line
# calling namd() functions (or may be some other unnecessary calculations)

Nstates = len(params["states"])  # Total number of basis states
inp_dir = os.getcwd()+"/out"     # this is the directory containing the input for this stage
                                 # it is the directory where pyxaid_core.namd() has written all
                                 # it output (raw output)
opt = 12                         # Defines the type of the averaging we want to do. Possible values:
                                 # 1 - average over intial conditions, independnetly for each state
                                 # 2 - sum the averages for groups of states (calculations with opt=1 must
                                 # already be done). One can do this for different groups of states without 
                                 # recomputing initial-conditions averages - they stay the same
                                 # 12 - do the steps 1 and 2 one after another

# Define the groups of states for which we want to know the total population as a function of time
MS = []
#for i in range(0,Nstates):
#    MS.append([i])   # In our case - each group of states (macrostate) contains only a single basis configuration
                     # (microstate)
MS.append([0])
MS.append([1])


res_dir = os.getcwd()+"/macro"  # Hey! : you need to create this folder in the current directory
                                # This is where the averaged results will be written

# Finally, run the averaging
average.average(params["namdtime"],Nstates,params["iconds"],opt,MS,inp_dir,res_dir)


