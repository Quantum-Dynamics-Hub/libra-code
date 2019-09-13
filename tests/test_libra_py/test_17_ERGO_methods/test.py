import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from libra_py import ERGO_methods
from libra_py import units
from libra_py import data_read



S = ERGO_methods.get_mtx_matrices("S_matrix.mtx")
print("S"); S.show_matrix()

case = 1  # 0 - non-polarized, 1 - polarized


if case==0:
    [e_a], [nocc, nvirt] = ERGO_methods.read_spectrum_restricted() 
    print(e_a, nocc, nvirt)

    [mo_a] = ERGO_methods.read_mo_restricted(nocc, nvirt) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()

    [mo_a] = ERGO_methods.read_mo_restricted(nocc, nvirt, [0, 1]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()

    [mo_a] = ERGO_methods.read_mo_restricted(nocc, nvirt, [-1, 0, 1]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()

    [mo_a] = ERGO_methods.read_mo_restricted(nocc, nvirt, [-1, 2]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()

elif case==1:
    [e_a, e_b], [nocc, nvirt] = ERGO_methods.read_spectrum_unrestricted() 
    print(e_a, e_b, nocc, nvirt)

    [mo_a, mo_b] = ERGO_methods.read_mo_unrestricted(nocc, nvirt) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()
    print("mo_b.H() * S * mo_b = ");  (mo_b.H() * S * mo_b).show_matrix()

    [mo_a, mo_b] = ERGO_methods.read_mo_unrestricted(nocc, nvirt, [0, 1]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()
    print("mo_b.H() * S * mo_b = ");  (mo_b.H() * S * mo_b).show_matrix()

    [mo_a, mo_b] = ERGO_methods.read_mo_unrestricted(nocc, nvirt, [-1, 0, 1]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()
    print("mo_b.H() * S * mo_b = ");  (mo_b.H() * S * mo_b).show_matrix()

    [mo_a, mo_b] = ERGO_methods.read_mo_unrestricted(nocc, nvirt, [-1, 2]) 
    print("mo_a.H() * S * mo_a = ");  (mo_a.H() * S * mo_a).show_matrix()
    print("mo_b.H() * S * mo_b = ");  (mo_b.H() * S * mo_b).show_matrix()


