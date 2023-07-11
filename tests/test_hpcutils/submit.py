"""
 This example will "submit" several simple jobs
"""

import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import hpc_utils


for i in range(5):
    os.system("mkdir res%i" %i)
    os.system("cp run.py res%i" %i)
    hpc_utils.substitute("run.py", "res%i/run.py" % (i), { "var=": "var = %i" % (i) })

    os.chdir("res%i" % i)
    os.system("python run.py")
    os.chdir("../")

