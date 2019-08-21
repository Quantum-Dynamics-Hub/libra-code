import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
from libra_py import ERGO_methods
from libra_py import units
from libra_py.workflows.nbra import step2_ergoscf


def run(EXE, COORDS):
    inp = """#!bin/sh
%s << EOINPUT > /dev/null
spin_polarization = 0
molecule_inline Angstrom
%sEOF
basis = "STO-3G"
use_simple_starting_guess=1
scf.create_mtx_files_F = 1
scf.create_mtx_file_S = 1
XC.sparse_mode = 1
run "LDA"
EOINPUT
""" % (EXE, COORDS)
    return inp

def run1(EXE, COORDS):
    inp = """#!bin/sh
%s << EOINPUT > /dev/null
spin_polarization = 0
molecule_inline Angstrom
%sEOF
basis = "STO-3G"
use_simple_starting_guess=1
scf.create_mtx_file_S = 1
scf.create_mtx_files_S_and_quit = 1
XC.sparse_mode = 1
run "LDA"
EOINPUT
""" % (EXE, COORDS)
    return inp


params = {"EXE":"/mnt/c/cygwin/home/Alexey-user/Soft/ergo-3.8/source/ergo", "md_file":"../step1/md.xyz",
          "isnap":0, "fsnap":5000, "dt": 1.0 * units.fs2au
         }

os.system("mkdir res")
step2_ergoscf.run_step2(params, run, run1)
