import pytest

from liblibra_core import *

import numpy as np
import libra_py.data_conv as data_conv

# Sets for testing edc_rates
set1 = [
[ np.array([[1.0 + 0.0j, 0.5 + 0.3j], [0.5 - 0.3j, 2.5 + 0.0j]]), 1.0, 1.0, 0.1, 0, np.array([[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]]) ],
[ np.array([[2.0 + 0.0j, 1.0 + 0.2j], [1.0 - 0.2j, 3.0 + 0.0j]]), 1.0, 1.0, 0.1, 0, np.array([[0.0, 0.9090909090909091], [0.9090909090909091, 0.0]]) ],
[ np.array([[0.5 + 0.0j, 0.3 + 0.1j, 0.2 - 0.1j], [0.3 - 0.1j, 1.0 + 0.0j, 0.4 + 0.0j], [0.2 + 0.1j, 0.4 + 0.0j, 1.5 + 0.0j]]), 1.0, 1.0, 0.1, 0,  
  np.array([[0.0, 0.45454545454545453, 0.9090909090909091], [0.45454545454545453 , 0.0, 0.45454545454545453], [ 0.9090909090909091, 0.45454545454545453, 0.0]]) ]
]

set2 = [
[ np.array([[0.8 + 0.1j], [0.6 - 0.2j]]), np.array([[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]]), np.array([[1.8333333333333335], [1.128205128205128]]) ],
[ np.array([[0.8 + 0.1j], [0.6 - 0.2j]]), np.array([[0.0, 0.9090909090909091], [0.9090909090909091, 0.0]]), np.array([[2.7500000000000000], [1.6923076923076918]]) ],
[ np.array([[1/np.sqrt(2)], [1/np.sqrt(2)]]), np.array([[0.0, 1.3636363636363636], [1.3636363636363636, 0.0]]), np.array([[1.466666666666667], [1.466666666666667]]) ],
[ np.array([[1/np.sqrt(2)], [1/np.sqrt(2)]]), np.array([[0.0, 0.9090909090909091], [0.9090909090909091, 0.0]]), np.array([[2.2000000000000006], [2.2000000000000006]]) ]
]

class TestSpecialFunctions: 

    @pytest.mark.parametrize(('Hvib', 'Ekin', 'C_param', 'eps_param', 'isNBRA', 'decoh_rates_expt'), set1)
    def test_1(self, Hvib, Ekin, C_param, eps_param, isNBRA, decoh_rates_expt):
        nst = Hvib.shape[0]
        Hvib = data_conv.nparray2CMATRIX(Hvib)

        decoh_rates = edc_rates(Hvib, Ekin, C_param, eps_param, isNBRA)
        decoh_rates = data_conv.MATRIX2nparray(decoh_rates)

        for ist in range(nst):
            for jst in range(nst):
                assert abs(decoh_rates[ist, jst] - decoh_rates_expt[ist, jst]) < 1e-10

    @pytest.mark.parametrize(('Coeff', 'rates', 'tau_m_expt'), set2)
    def test_2(self, Coeff, rates, tau_m_expt):
        nst = Coeff.shape[0]
        
        Coeff = data_conv.nparray2CMATRIX(Coeff)
        rates = data_conv.nparray2MATRIX(rates)

        tau_m = coherence_intervals(Coeff, rates)
        tau_m = data_conv.MATRIX2nparray(tau_m)

        for ist in range(nst):
          assert abs(tau_m[ist, 0] - tau_m_expt[ist, 0]) < 1e-10

