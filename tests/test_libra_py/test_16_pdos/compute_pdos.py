import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import pdos as pdos


#Example of usage:
spd = ["s","p","d"]
PT = ["Cs", "Sn", "Br"]
E_f = 6.5084  # Fermi energy

"""
Cs = [spd, range(1,25), PT]
Sn = [spd, range(25,31), PT]
Br = [spd, range(31,67), PT]
projections = [ Cs, Sn, Br ]
"""

"""
Cs = [spd, range(1,100), ["Cs"] ]
Sn = [spd, range(1,100), ["Sn"] ]
Br = [spd, range(1,100), ["Br"] ]
projections = [ Cs, Sn, Br ]
"""

#pdos.QE_pdos("pdos/Cs4SnBr6_T100.pdos.pdos_atm#", -10.0, 10.0, 0.1, projections, E_f, "pdos.txt", 1, 0.01, 0.1)



Sns = [["s"], range(1,100), ["Sn"] ]
Snp = [["p"], range(1,100), ["Sn"] ]
Snd = [["d"], range(1,100), ["Sn"] ]
projections = [ Sns, Snp, Snd ]

pdos.QE_pdos("pdos/Cs4SnBr6_T100.pdos.pdos_atm#", -10.0, 10.0, 0.1, projections, E_f, "pdos.txt", 1, 0.01, 0.1)



