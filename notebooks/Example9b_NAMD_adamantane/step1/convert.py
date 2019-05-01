import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import DFTB_methods


DFTB_methods.xyz_traj2gen_sp("adamantane.xyz", "adamantane.gen", 0, "C")
