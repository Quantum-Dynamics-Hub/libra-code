import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
    

from libra_py import hpc_utils
from libra_py import data_read
from libra_py import data_outs
from libra_py import units
from libra_py import QE_methods
from libra_py.workflows.nbra import step2


print( os.getcwd())



PP_dir = os.getcwd()+"/PP/"

scf_in = """&CONTROL
  calculation = 'scf',
  dt = 20.67055,
  nstep = 10,
  pseudo_dir = '%s',
  outdir = './',
  prefix = 'x0',
  disk_io = 'low',
  wf_collect = .true.
/

&SYSTEM
  ibrav = 0,
  celldm(1) = 1.89,
  nat = 4,
  ntyp = 2,
  nspin = 2,
  starting_magnetization(1) = 0.1,
  nbnd = 40,
  ecutwfc = 40,
  tot_charge = 0.0,
  occupations = 'smearing',
  smearing = 'gaussian',
  degauss = 0.005,
  nosym = .true.,
/

&ELECTRONS
  electron_maxstep = 300,
  conv_thr = 1.D-5,
  mixing_beta = 0.45,
/

&IONS
  ion_dynamics = 'verlet',
  ion_temperature = 'andersen',
  tempw = 300.00 ,
  nraise = 1,
/


ATOMIC_SPECIES
 Cd 121.411  Cd.pbe-n-rrkjus_psl.1.0.0.UPF
 Se 78.96    Se.pbe-dn-rrkjus_psl.1.0.0.UPF


K_POINTS automatic
 1 1 1 0 0 0

CELL_PARAMETERS (alat=  1.89000000)
   4.716986504  -0.015512615  -0.002400656
  -2.371926710   4.062829845  -0.000273730
  -0.002552594  -0.001387965   8.436361230

""" % (PP_dir)

exp_in = """&inputpp
  prefix = 'x0',
  outdir = './',
  pseudo_dir = '%s',
  psfile(1) = 'Cd.pbe-n-rrkjus_psl.1.0.0.UPF',
  psfile(2) = 'Se.pbe-dn-rrkjus_psl.1.0.0.UPF',
  single_file = .FALSE.,
  ascii = .TRUE.,
  uspp_spsi = .FALSE.,
/

""" % (PP_dir)

f = open("x0.scf.in", "w")
f.write(scf_in)
f.close()

f = open("x0.exp.in", "w")
f.write(exp_in)
f.close()


#print scf_in
#print exp_in



# Remove the previous results and temporary working directory from the previous runs
os.system("rm -r res")
os.system("rm -r wd")

# Create the new results directory
os.system("mkdir res")
rd = os.getcwd()+"/res"          # where all the results stuff will go




submit_str = """#!/bin/sh
#SBATCH --partition=valhalla --qos=valhalla
#SBATCH --clusters=chemistry
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=5000
#SBATCH --mail-user=alexeyak@buffalo.edu

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory="$SLURM_SUBMIT_DIR

NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS

module load intel/16.0
module load intel-mpi/5.1.1
module load mkl/11.3
module load espresso


#The PMI library is necessary for srun
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# These will be assigned automatically, leave them as they are
param1=
param2=

# This is invocation of the scripts which will further handle NA-MD calclculations on the NAC calculation step
# NOTE: minband - starting from 1
#       maxband - is included

python -c \"from libra_py.workflows.nbra import step2
params = {}
params[\\"EXE\\"] = \\"pw.x\\" 
params[\\"EXE_EXPORT\\"] = \\"pw_export.x\\"
params[\\"BATCH_SYSTEM\\"] = \\"srun\\"
params[\\"NP\\"] = 4
params[\\"start_indx\\"] = $param1
params[\\"stop_indx\\"] = $param2
params[\\"dt\\"] = %8.5f
params[\\"prefix0\\"] = \\"x0.scf\\" 
params[\\"nac_method\\"] = 1
params[\\"minband\\"] = 20
params[\\"maxband\\"] = 39
params[\\"minband_soc\\"] = 20
params[\\"maxband_soc\\"] = 39
params[\\"compute_Hprime\\"] = True
params[\\"wd\\"] = \\"wd\\"
params[\\"rd\\"] = \\"%s\\"
params[\\"verbosity\\"] = 0
step2.run(params)
\"

""" % ( 1.0*units.fs2au, os.getcwd()+"/res" )


f = open("submit_templ.slm", "w")
f.write(submit_str)
f.close()

#print submit_str


iinit = 0
ifinal = 10

QE_methods.out2inp("x0.md.out","x0.scf.in","wd","x0.scf", iinit, ifinal,1)

os.system("cp submit_templ.slm wd")
os.system("cp x0.exp.in wd")


os.chdir("wd")

tot_nsteps = 10
nsteps_per_job = 5
hpc_utils.distribute(0,tot_nsteps,nsteps_per_job,"submit_templ.slm",["x0.exp.in"],["x0.scf"],2)

os.chdir("../")


print( os.getcwd())

