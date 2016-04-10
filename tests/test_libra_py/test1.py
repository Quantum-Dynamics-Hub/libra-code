import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

label = ["H"]
R = [VECTOR()]
g = [VECTOR()]
rnd = Random()
T = 300.0
sigma = 0.1
df = 1

x = init_system.init_system(label, R, g, rnd, T, sigma, df)

