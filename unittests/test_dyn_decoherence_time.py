import pytest

import math
from liblibra_core import *

# Sets for testing edc_rates
set1 = [
[ [[1.0 + 0.0j, 0.5 + 0.3j], [0.5 - 0.3j, 1.5 + 0.0j]], 1.0e+20, 1.0, 0.1, 0, [[0.0, 0.5], [0.5, 0.0]] ], # When Ekin -> \infty, rates = |\Delta E| / \hbar
[ [[1.0 + 0.0j, 0.5 + 0.3j], [0.5 - 0.3j, 2.5 + 0.0j]], 1.0, 1.0, 0.1, 0, [[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]] ],
[ [[2.0 + 0.0j, 1.0 + 0.2j], [1.0 - 0.2j, 3.0 + 0.0j]], 1.0, 1.0, 0.1, 0, [[0.0, 0.9090909090909091], [0.9090909090909091, 0.0]] ],
[ [[0.5 + 0.0j, 0.3 + 0.1j, 0.2 - 0.1j], [0.3 - 0.1j, 1.0 + 0.0j, 0.4 + 0.0j], [0.2 + 0.1j, 0.4 + 0.0j, 1.5 + 0.0j]], 1.0, 1.0, 0.1, 0,  
  [[0.0, 0.45454545454545453, 0.9090909090909091], [0.45454545454545453 , 0.0, 0.45454545454545453], [ 0.9090909090909091, 0.45454545454545453, 0.0]] ]
]

# Sets for testing coherence_intervals
set2 = [
[ [[1/math.sqrt(2)], [1/math.sqrt(2)]], [[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]], [[1.466666666666667], [1.466666666666667]] ], # equal populations
[ [[0.8 + 0.1j], [0.6 - 0.2j]], [[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]], [[1.8333333333333335], [1.128205128205128]] ],
[ [[0.8 + 0.1j], [0.6 - 0.2j]], [[0.0, 0.9090909090909091], [0.9090909090909091, 0.0]], [[2.7500000000000000], [1.6923076923076918]] ],
]

class TestSpecialFunctions: 

    @pytest.mark.parametrize(('_Hvib', '_Ekin', '_C_param', '_eps_param', '_isNBRA', 'res_expt'), set1)
    def test_1(self, _Hvib, _Ekin, _C_param, _eps_param, _isNBRA, res_expt):
        nst = len(_Hvib)

        Hvib = CMATRIX(nst, nst)
        for i in range(nst):
            for j in range(nst):
                Hvib.set(i, j, Py2Cpp_complex([_Hvib[i][j]])[0] )
        
        Ekin = Py2Cpp_double([_Ekin])[0]
        C_param = Py2Cpp_double([_C_param])[0]
        eps_param = Py2Cpp_double([_eps_param])[0]
        isNBRA = Py2Cpp_int([_isNBRA])[0]

        decoh_rates = edc_rates(Hvib, Ekin, C_param, eps_param, isNBRA)

        for i in range(nst):
            for j in range(nst):
                assert abs(decoh_rates.get(i, j) - res_expt[i][j]) < 1e-10

    @pytest.mark.parametrize(('_Coeff', '_rates', 'res_expt'), set2)
    def test_2(self, _Coeff, _rates, res_expt):
        nst = len(_Coeff)
        
        Coeff = CMATRIX(nst, 1)
        for i in range(nst):
            Coeff.set(i, 0, Py2Cpp_complex([_Coeff[i][0]])[0] )
        
        rates = MATRIX(nst, nst)
        for i in range(nst):
            for j in range(nst):
                rates.set(i, j, Py2Cpp_double([_rates[i][j]])[0] )

        tau_m = coherence_intervals(Coeff, rates)

        for i in range(nst):
          assert abs(tau_m.get(i, 0) - res_expt[i][0]) < 1e-10

