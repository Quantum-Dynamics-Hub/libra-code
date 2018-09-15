import math
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


def run_test():

    # Parameters
    inv_cm2ev = (1.0/8065.54468111324)
    ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
    inv_cm2Ha = inv_cm2ev * ev2Ha
    fs2au = (1.0/0.02419)   # 40 a.u. is 1 fs 

    # Test case: 3 frequences
    data = []
    dt = 1.0
    dw = 1.0
    wspan = 2000.0
    w1 = 500.0 * inv_cm2Ha
    w2 = 1400.0 * inv_cm2Ha
    w3 = 850.0 * inv_cm2Ha

    for it in xrange(1000):
        t = it * dt * fs2au
        data.append( VECTOR(math.sin(w1*t), math.cos(w2*t), math.sin(w3*t)) )
    
    acf_vector.recipe1(data, dt, wspan, dw)

run_test()
