#*********************************************************************************                     
#* Copyright (C) 2018-2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: models_LVC
   :platform: Unix, Windows
   :synopsis: This module implements Linear Vibronic Coupling (LVC) Hamiltonian
.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units


class tmp:
    pass    


def LVC(q, params, full_id=None):
    """
    Linear Vibronic Coupling Hamiltonian, 2-level, N-dim. problem

    Args:
        q ( MATRIX(ndof, 1) ): coordinates of the classical particles, ndof is an 
            arbitrary number of degrees of freedom (e.g. 3N, where N is the number of particles)
        params ( dictionary ): model parameters, should contain:

            * **params["Delta1"]** ( double ): energy minimum of the lower state [ units: Ha ]
            * **params["Delta2"]** ( double ): energy minimum of the upper state [ units: Ha ]
            * **params["omega"]** ( list on ndof doubles ): normal modes frequencies
                same for both electronic states [ units: Ha ]
            * **params["d1"]** ( list on ndof doubles ): electron-phonon couplings for 
                the lower state [ units: Ha/Bohr ]
            * **params["d2"]** ( list on ndof doubles ): electron-phonon couplings for 
                the upper state [ units: Ha/Bohr ]
            * **params["coup"]** ( list on ndof doubles ): electron-phonon couplings [ units: Ha/Bohr ]
            * **params["mass"]** ( list on ndof doubles ): masses of the normal modes [ units: amu ]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of ndof CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]


    Note: the Hamiltonian defined in

    Izmaylov, A. F.; Mendive-Tapia, D.; Bearpark, M. J.; Robb, M. A.; 
    Tully, J. C.; Frisch, M. J. JCP, 2011, 135, 234106

    uses the mass-transformed coordinates and momenta. 
    In the present formulation, the regular coordinates (not mass-transformed) are used.
    
    """

    # Define potential specific constants
    critical_params = [ "omega", "d1", "d2", "coup", "mass", "Delta1", "Delta2" ] 
    default_params = { }
    comn.check_input(params, default_params, critical_params)

    w = params["omega"]
    d1 = params["d1"]
    d2 = params["d2"]
    c = params["coup"]
    m = params["mass"]



    ndof = q.num_of_rows  # the number of nuclear DOFs

    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(0,ndof):
        obj.d1ham_dia.append( CMATRIX(2,2) )
        obj.dc1_dia.append( CMATRIX(2,2) )


    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]


    #=========== Energies & Derivatives ===============
    obj.ham_dia.set(0,0, params["Delta1"]*(1.0+0.0j))
    obj.ham_dia.set(1,1, params["Delta2"]*(1.0+0.0j))

    for n in range(0,ndof):
        # Diagonal bath
        q_n = q.get(n, indx)

        x = 0.5 * m[n] * w[n] * w[n]* q_n * q_n * (1.0+0.0j) 
        obj.ham_dia.add(0,0, x)
        obj.ham_dia.add(1,1, x)

        x = m[n] * w[n] * w[n] * q_n * (1.0+0.0j) 
        obj.d1ham_dia[n].add(0, 0, x)
        obj.d1ham_dia[n].add(1, 1, x)


        # Diagonal el-ph coupling
        x = math.sqrt(m[n]) * q_n * (1.0+0.0j) 
        obj.ham_dia.add(0,0, x*d1[n])
        obj.ham_dia.add(1,1, x*d2[n])

        x = math.sqrt(m[n])*(1.0+0.0j) 
        obj.d1ham_dia[n].add(0, 0, x*d1[n])
        obj.d1ham_dia[n].add(0, 0, x*d2[n])


        # Off-diagonal
        x = math.sqrt(m[n]) * c[n] * (1.0+0.0j)
        obj.ham_dia.add(0,1, x*q.get(n,0) )
        obj.ham_dia.add(1,0, x*q.get(n,0) )

        obj.d1ham_dia[n].add(0, 0, x)
        obj.d1ham_dia[n].add(0, 0, x)


    return obj



def get_LVC_set1():
    """

    Parameters for **Fulvene** molecule
    References: 

        1. Izmaylov, A. F.; Mendive-Tapia, D.; Bearpark, M. J.; Robb, M. A.; 
            Tully, J. C.; Frisch, M. J. JCP, 2011, 135, 234106

        2. Sun, X.; Geva, E. J. Chem. Phys. 2016, 144, 244105 (files)

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

            * **params["omega_DA"]** ( double ): donor-acceptor energy gap [ units: Ha ] 
                TODO: how is this different from Delta2 - Delta1?
            * **params["Delta1"]** ( double ): energy minimum of the lower state [ units: Ha ]
            * **params["Delta2"]** ( double ): energy minimum of the upper state [ units: Ha ]
            * **params["Er"]** ( double ): reorganization energy [ units: Ha ]
            * **params["omega"]** ( list on ndof doubles ): normal modes frequencies
                same for both electronic states [ units: Ha ]
            * **params["d1"]** ( list on ndof doubles ): electron-phonon couplings for 
                the lower state [ units: Ha/Bohr ]
            * **params["d2"]** ( list on ndof doubles ): electron-phonon couplings for 
                the upper state [ units: Ha/Bohr ]
            * **params["coup"]** ( list on ndof doubles ): electron-phonon couplings [ units: Ha/Bohr ]

    """

    params = {}
    params["omega_DA"]  = 0.0989  # Donor-Acceptor energy gap, a.u.
    params["Er"] = 0.0887         # Reorganization energy, a.u.
    params["Delta1"] = 0.0
    params["Delta2"] = 1.8764e-01

    params["omega"] = [9.8662e-04, 1.7237e-03, 2.4190e-03, 2.8634e-03, 3.3546e-03,
                       3.4833e-03, 3.6851e-03, 3.8835e-03, 4.3459e-03, 4.1842e-03,
                       4.3828e-03, 4.4617e-03, 4.6002e-03, 5.1482e-03, 5.1310e-03,
                       5.7261e-03, 5.7567e-03, 6.7668e-03, 7.7550e-03, 7.2129e-03,
                       7.6025e-03, 7.9270e-03, 8.2881e-03, 8.5797e-03, 1.6669e-02,
                       1.7018e-02, 1.7046e-02, 1.7134e-02, 1.7176e-02, 1.7397e-02 ]

    params["d1"] =    [0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0        ]

    params["d2"] =    [0.0       ,-9.6755e-09, 0.0       , 0.0       , 6.4918e-03,
                       0.0       , 0.0       , 0.0       , 2.9431e-07, 0.0       ,
                       0.0       , 0.0       , 7.1471e-03, 3.9572e-08, 6.4988e-03,
                       1.9246e-06,-7.4781e-03,-6.4437e-07, 1.1289e-07, 7.2115e-03,
                       5.5555e-03,-1.8578e-02, 6.9979e-07, 2.6226e-02,-3.6746e-03,
                      -4.7285e-07, 7.4650e-04,-7.3343e-07,-6.4173e-04, 1.1643e-07 ]

    params["coup"] =  [0.0       ,-2.2284e-04, 0.0       , 0.0       ,-3.6072e-09,
                       0.0       , 0.0       , 0.0       , 6.5332e-03, 0.0       ,
                       0.0       , 0.0       ,-1.5055e-07, 5.2140e-03,-3.0083e-08,
                       6.0433e-04, 1.2886e-07, 6.2708e-03, 1.0772e-02, 3.7441e-07,
                       1.4169e-07,-2.3477e-08,-3.8450e-03,-4.2660e-08, 1.5939e-09,
                       5.4285e-05,-2.5890e-08,-3.9978e-04, 7.9818e-08, 2.2991e-05 ]

    return params



def get_LVC_set1b():
    """

    Parameters for **Fulvene** molecule
    From the NEFGRL_Fulvene.cpp code

    References: 

        1. Sun, X.; Geva, E. J. Chem. Phys. 2016, 144, 244105 (code!)

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

            * **params["omega_DA"]** ( double ): donor-acceptor energy gap [ units: Ha ] 
                TODO: how is this different from Delta2 - Delta1?
            * **params["Delta1"]** ( double ): energy minimum of the lower state [ units: Ha ]
            * **params["Delta2"]** ( double ): energy minimum of the upper state [ units: Ha ]
            * **params["Er"]** ( double ): reorganization energy [ units: Ha ]
            * **params["omega"]** ( list on ndof doubles ): normal modes frequencies
                same for both electronic states [ units: Ha ]
            * **params["d1"]** ( list on ndof doubles ): electron-phonon couplings for 
                the lower state [ units: Ha/Bohr ]
            * **params["d2"]** ( list on ndof doubles ): electron-phonon couplings for 
                the upper state [ units: Ha/Bohr ]
            * **params["coup"]** ( list on ndof doubles ): electron-phonon couplings [ units: Ha/Bohr ]



    """

    params = {}
    params["omega_DA"]  = 0.0989  # Donor-Acceptor energy gap, a.u.
    params["Er"] = 0.0887         # Reorganization energy, a.u.
    params["Delta1"] = 0.0
    params["Delta2"] = 1.8764e-01

    params["omega"] = [0.00098662,  0.0017237,   0.002419,  0.0028634,  0.0033546,
                        0.0034833,  0.0036851,  0.0038835,  0.0043459,  0.0041842,
                        0.0043828,  0.0044617,  0.0046002,  0.0051482,   0.005131,
                        0.0057261,  0.0057567,  0.0067668,   0.007755,  0.0072129,
                        0.0076025,   0.007927,  0.0082881,  0.0085797,   0.016669,
                         0.017018,   0.017046,   0.017134,   0.017176,   0.017397 ]

    params["d1"] =    [0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0        ]

    params["d2"] =    [       0.0, -4.017e-10,        0.0,        0.0,   0.000376,
                              0.0,        0.0,        0.0, 1.9402e-08,        0.0,
                              0.0,        0.0, 0.00048475, 2.8393e-09, 0.00046551,
                       1.4564e-07,-0.00056739,-5.3006e-08, 9.9417e-09, 0.00061246,
                        0.0004844, -0.0016541, 6.3709e-08,  0.0024292,-0.00047442,
                      -6.1685e-08, 9.7464e-05,-9.6004e-08,-8.4104e-05, 1.5356e-08 ]

    params["coup"] =  [       0.0,-9.2519e-06,        0.0,        0.0,-2.0893e-10,
                              0.0,        0.0,        0.0, 0.00043069,        0.0,
                              0.0,        0.0,-1.0211e-08, 0.00037411,-2.1549e-09,
                        4.573e-05, 9.7768e-09, 0.00051583,  0.0009486, 3.1798e-08,
                       1.2354e-08,-2.0903e-09,-0.00035004,-3.9514e-09, 2.0579e-10,
                       7.0815e-06,-3.3803e-09, -5.233e-05, 1.0461e-08, 3.0324e-06 ]



    return params
   

   

def get_LVC_set2():
    """

    Parameters for **2,6-bis(methylene)adamantyl (BMA)** radical cation
    References: 

        1. Izmaylov, A. F.; Mendive-Tapia, D.; Bearpark, M. J.; Robb, M. A.; 
            Tully, J. C.; Frisch, M. J. JCP, 2011, 135, 234106

        2. Sun, X.; Geva, E. J. Chem. Phys. 2016, 144, 244105 (files)

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

            * **params["omega_DA"]** ( double ): donor-acceptor energy gap [ units: Ha ] 
                TODO: how is this different from Delta2 - Delta1?
            * **params["Delta1"]** ( double ): energy minimum of the lower state [ units: Ha ]
            * **params["Delta2"]** ( double ): energy minimum of the upper state [ units: Ha ]
            * **params["Er"]** ( double ): reorganization energy [ units: Ha ]
            * **params["omega"]** ( list on ndof doubles ): normal modes frequencies
                same for both electronic states [ units: Ha ]
            * **params["d1"]** ( list on ndof doubles ): electron-phonon couplings for 
                the lower state [ units: Ha/Bohr ]
            * **params["d2"]** ( list on ndof doubles ): electron-phonon couplings for 
                the upper state [ units: Ha/Bohr ]
            * **params["coup"]** ( list on ndof doubles ): electron-phonon couplings [ units: Ha/Bohr ]


    """

    params = {}
    params["omega_DA"]  = 0.0     # Donor-Acceptor energy gap, a.u.
    params["Er"] = 0.0297         # Reorganization energy, a.u.
    params["Delta1"] = 0.0
    params["Delta2"] = 2.9255e-02

    params["omega"] =[ 1.4364e-04, 7.0042e-04, 1.3609e-03, 1.4446e-03, 1.5662e-03,
                       1.7429e-03, 1.9084e-03, 1.9872e-03, 2.0680e-03, 2.1338e-03,
                       2.2608e-03, 2.9137e-03, 3.0567e-03, 3.1336e-03, 3.5698e-03,
                       3.5937e-03, 3.6504e-03, 3.7280e-03, 3.9665e-03, 4.2489e-03,
                       4.3519e-03, 4.4641e-03, 4.5020e-03, 4.5767e-03, 4.8162e-03,
                       4.8881e-03, 4.9876e-03, 5.0317e-03, 5.0655e-03, 5.3204e-03,
                       5.3536e-03, 5.5110e-03, 5.5488e-03, 5.7590e-03, 5.8294e-03,
                       5.9130e-03, 6.0599e-03, 6.1243e-03, 6.2667e-03, 6.4994e-03,
                       6.7154e-03, 6.7344e-03, 6.8892e-03, 6.9729e-03, 6.9762e-03,
                       7.1230e-03, 7.1467e-03, 7.2810e-03, 7.3408e-03, 7.4179e-03,
                       7.4193e-03, 7.4574e-03, 7.4792e-03, 7.5643e-03, 7.5861e-03,
                       7.8827e-03, 7.9413e-03, 8.1991e-03, 8.2369e-03, 8.2473e-03,
                       8.2744e-03, 8.7380e-03, 1.6495e-02, 1.6496e-02, 1.6500e-02,
                       1.6500e-02, 1.6556e-02, 1.6590e-02, 1.6717e-02, 1.6724e-02,
                       1.6779e-02, 1.6789e-02, 1.7025e-02, 1.7032e-02, 1.7032e-02,
                       1.7038e-02, 1.7316e-02, 1.7437e-02]

    params["d1"] =   [ 0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       ]

    params["d2"] =   [ 9.5316e-08,-6.1610e-09, 2.1465e-08, 1.3506e-08,-1.2410e-08,
                      -9.7421e-05, 7.4580e-08, 5.1857e-08, 1.9686e-07,-7.8621e-05,
                       2.6906e-08, 6.8739e-04, 2.5550e-08, 1.9812e-07,-8.6642e-05,
                      -5.5297e-09,-5.6090e-09, 1.5535e-07,-7.0806e-04, 7.5569e-09,
                      -1.6777e-09, 6.9802e-08, 1.1450e-07,-1.4488e-03, 1.1084e-08,
                       2.2397e-08, 1.0841e-07,-2.4316e-03, 2.3064e-07, 1.2244e-04,
                      -5.4065e-09, 1.8170e-08, 1.8416e-08, 8.6517e-08, 4.6066e-08,
                       4.2798e-04, 4.3259e-08, 3.1798e-09,-3.0551e-08,-2.5659e-03,
                       1.0659e-02,-5.0572e-08, 6.6957e-08, 3.9842e-08,-2.4300e-08,
                       3.0742e-03, 3.4232e-08,-5.2202e-08,-8.4530e-08, 3.9075e-08,
                       2.3079e-08, 5.7213e-03, 1.5880e-07,-4.9429e-07,-9.5236e-03,
                      -4.3176e-03, 5.2027e-03,-9.5068e-09,-5.1084e-08, 3.1476e-08,
                       2.1445e-04, 1.1784e-02,-5.5719e-08,-1.7616e-07, 1.0770e-07,
                       2.2482e-05, 9.4927e-04,-7.3844e-04,-1.8251e-07, 1.4299e-04,
                       1.4446e-09, 1.0025e-04, 2.9834e-08, 2.1950e-07,-4.7783e-08,
                       9.4639e-05,-8.9160e-08, 1.3428e-10]

    params["coup"] = [-2.1202e-09,-6.5125e-10,-1.0931e-05, 1.3263e-08,-3.4265e-10,
                      -9.3686e-10,-9.2689e-09,-3.0003e-04, 2.8083e-09,-7.5714e-09,
                       1.9786e-04,-5.6040e-09,-1.9112e-09,-8.9604e-10,-9.0550e-10,
                       2.6481e-04, 4.7121e-09, 1.3096e-09,-4.7597e-09, 1.1808e-08,
                       1.6810e-09, 2.2826e-09,-1.7983e-09, 2.8143e-09, 4.1939e-06,
                      -2.9435e-09, 2.2324e-04,-2.9482e-09, 1.2870e-09,-4.5466e-10,
                      -7.5413e-09, 2.9948e-09, 2.4552e-10,-1.0830e-09,-1.4443e-05,
                       1.9933e-10,-2.4312e-09, 3.1468e-09,-2.6851e-04, 1.0466e-09,
                      -2.4879e-09,-1.9627e-04, 1.2039e-11,-1.9576e-09, 3.0448e-09,
                      -4.8520e-10,-5.7719e-09, 2.7622e-04, 6.2422e-10, 7.0792e-04,
                      -3.4547e-08, 3.2337e-09, 4.1860e-09,-4.4150e-09,-1.4466e-09,
                       1.2159e-09,-3.4545e-10,-2.3328e-05, 6.0943e-09, 2.9680e-09,
                       2.3591e-09, 1.6218e-10, 2.2328e-04,-2.7163e-07, 3.7086e-07,
                      -1.5156e-07,-3.0636e-10, 4.7115e-11, 2.6172e-10, 3.4384e-10,
                      -1.4759e-10,-6.3461e-10, 2.1464e-05,-1.0486e-08, 3.4116e-08,
                       8.3360e-09,-2.1252e-11, 1.9538e-10]

    return params
   

def get_LVC_set3():
    """

    Parameters for **2-methylene-6-isopropylidene-adamantyl (MIA)** radical cation
    References: 

        1. Izmaylov, A. F.; Mendive-Tapia, D.; Bearpark, M. J.; Robb, M. A.; 
            Tully, J. C.; Frisch, M. J. JCP, 2011, 135, 234106

        2. Sun, X.; Geva, E. J. Chem. Phys. 2016, 144, 244105 (files)

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

            * **params["omega_DA"]** ( double ): donor-acceptor energy gap [ units: Ha ] 
                TODO: how is this different from Delta2 - Delta1?
            * **params["Delta1"]** ( double ): energy minimum of the lower state [ units: Ha ]
            * **params["Delta2"]** ( double ): energy minimum of the upper state [ units: Ha ]
            * **params["Er"]** ( double ): reorganization energy [ units: Ha ]
            * **params["omega"]** ( list on ndof doubles ): normal modes frequencies
                same for both electronic states [ units: Ha ]
            * **params["d1"]** ( list on ndof doubles ): electron-phonon couplings for 
                the lower state [ units: Ha/Bohr ]
            * **params["d2"]** ( list on ndof doubles ): electron-phonon couplings for 
                the upper state [ units: Ha/Bohr ]
            * **params["coup"]** ( list on ndof doubles ): electron-phonon couplings [ units: Ha/Bohr ]



    """

    params = {}
    params["omega_DA"]  = 0.0250  # Donor-Acceptor energy gap, a.u.
    params["Er"] = 0.0274         # Reorganization energy, a.u.
    params["Delta1"] = 0.0
    params["Delta2"] = 2.4615e-03

    params["omega"] = [3.2508e-05, 1.7800e-04, 3.3693e-04, 4.1833e-04, 6.8418e-04,
                       1.0397e-03, 1.3946e-03, 1.4123e-03, 1.4261e-03, 1.6660e-03,
                       1.8807e-03, 2.0590e-03, 2.0805e-03, 2.1475e-03, 2.2491e-03,
                       2.2329e-03, 2.4462e-03, 2.8244e-03, 3.0722e-03, 3.1405e-03,
                       3.4162e-03, 3.6092e-03, 3.7397e-03, 3.9021e-03, 4.3558e-03,
                       4.4143e-03, 4.4535e-03, 4.5541e-03, 4.7120e-03, 4.8109e-03,
                       4.8773e-03, 4.9855e-03, 5.0807e-03, 5.2753e-03, 5.3127e-03,
                       5.3129e-03, 5.3464e-03, 5.3726e-03, 5.5480e-03, 5.7395e-03,
                       5.7591e-03, 5.8360e-03, 5.9218e-03, 5.9608e-03, 6.0464e-03,
                       6.0585e-03, 6.2956e-03, 6.4666e-03, 6.5128e-03, 6.7295e-03,
                       6.7562e-03, 6.9176e-03, 6.9796e-03, 7.0121e-03, 7.1426e-03,
                       7.3223e-03, 7.3390e-03, 7.3420e-03, 7.3924e-03, 7.5859e-03,
                       7.4797e-03, 7.4872e-03, 7.5688e-03, 7.5837e-03, 7.8091e-03,
                       7.8638e-03, 7.9243e-03, 8.2062e-03, 8.2415e-03, 8.2545e-03,
                       8.2611e-03, 8.2762e-03, 8.2851e-03, 8.3019e-03, 8.3351e-03,
                       8.7153e-03, 1.6243e-02, 1.6246e-02, 1.6513e-02, 1.6494e-02,
                       1.6497e-02, 1.6499e-02, 1.6635e-02, 1.6791e-02, 1.6795e-02,
                       1.6812e-02, 1.6824e-02, 1.7022e-02, 1.7028e-02, 1.7030e-02,
                       1.7035e-02, 1.7042e-02, 1.7041e-02, 1.7132e-02, 1.7152e-02,
                       1.7487e-02]

    params["d1"] =   [ 0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0       , 0.0       , 0.0       , 0.0       , 0.0       ,
                       0.0 ]

    params["d2"] =   [-9.3059e-07, 8.4944e-06, 7.4335e-08, 2.0297e-06,-3.7710e-06,
                       2.1909e-08,-4.0926e-07, 2.9569e-04, 4.3562e-06, 8.9344e-08,
                      -3.5440e-07, 1.3791e-08, 5.6661e-06,-3.0083e-06,-3.2434e-07,
                       2.5619e-04, 1.4481e-08,-1.0063e-03, 1.8627e-07,-4.6159e-08,
                       3.4012e-04, 1.8648e-07,-6.1249e-08,-9.3485e-04,-1.8160e-07,
                      -8.7524e-04,-3.3247e-08, 7.0389e-09, 1.9783e-03, 1.6186e-06,
                      -1.1568e-07, 9.2399e-06,-1.2792e-03,-1.4293e-05,-9.4604e-07,
                       3.5716e-04,-5.7517e-08,-6.0132e-08, 2.4542e-08, 9.2559e-07,
                       3.8842e-07,-1.8082e-06,-1.9601e-04, 4.6979e-03,-5.0499e-08,
                       5.0949e-07, 8.6207e-07, 1.3872e-08,-4.7376e-03,-2.6562e-05,
                       8.3913e-03, 3.5207e-08,-1.0936e-07,-2.9351e-08, 3.8685e-08,
                      -3.0709e-05, 2.4427e-07,-7.3767e-05,-1.1101e-09, 1.5333e-05,
                      -2.1237e-06, 4.8574e-03, 1.0793e-06,-9.4477e-03, 7.8155e-09,
                       3.4019e-03, 6.3308e-03,-4.9945e-06,-1.3774e-08, 4.8773e-08,
                       1.1593e-07, 2.3460e-04, 4.0232e-04,-1.2460e-06,-1.5580e-06,
                       1.0918e-02,-1.5678e-06, 9.8301e-04,-1.9782e-05,-3.9410e-07,
                      -7.4772e-07, 7.9329e-05,-5.3783e-04, 6.9299e-07,-2.4104e-04,
                       2.0601e-07,-1.8099e-04,-2.1600e-06,-1.8077e-07, 1.1563e-07,
                      -2.0867e-05, 1.0952e-05, 3.6633e-08,-1.4219e-07,-4.8473e-04,
                      -6.9929e-08 ]

    params["coup"] = [-7.5605e-08,-2.6194e-04, 1.2815e-08,-8.9199e-09, 1.7312e-04,
                      -1.3217e-07, 2.2498e-08,-1.8086e-08, 5.4675e-06,-4.4629e-08,
                      -5.1885e-09,-1.9681e-09, 2.0204e-07, 3.0435e-04,-2.1092e-04,
                       4.6172e-07,-1.0543e-08,-1.1249e-06,-6.0323e-09, 1.3708e-08,
                       2.1967e-07,-3.1861e-09,-1.2269e-08,-1.1196e-06,-2.3180e-10,
                      -4.2942e-07,-8.8920e-09,-1.5255e-08, 1.5843e-06,-1.4916e-05,
                       6.8481e-10,-1.9256e-04,-1.2818e-06, 1.6108e-04, 5.6228e-08,
                       8.4848e-07, 4.6372e-09, 2.3664e-08, 2.2254e-08,-2.9472e-09,
                      -1.0357e-08, 7.5914e-05,-5.5114e-07, 5.6788e-06,-8.5787e-11,
                      -2.2638e-09,-2.5682e-04,-5.2163e-08,-5.4077e-06,-1.8016e-04,
                       8.5315e-06, 1.9762e-08, 1.3012e-09,-1.5839e-08, 1.0844e-09,
                       2.5769e-04, 8.2685e-09,-3.1629e-06, 1.4173e-07,-6.6751e-04,
                      -2.5883e-08, 6.3272e-06,-3.3997e-09,-1.0092e-05, 1.9688e-09,
                       3.8059e-06, 7.1273e-06, 3.0660e-05, 2.7060e-09,-1.5311e-09,
                      -1.1344e-08, 1.0562e-06, 1.5417e-06, 5.3452e-05, 2.1900e-09,
                       1.1678e-05,-3.2130e-08, 1.3009e-06,-2.2548e-04, 4.9117e-07,
                      -5.7472e-07,-2.7713e-05,-8.7741e-07, 9.1542e-10,-4.1811e-08,
                      -4.8796e-08, 2.2447e-07, 1.7331e-05,-1.4452e-08,-2.7777e-08,
                      -1.0289e-06,-6.2780e-05,-2.7487e-06, 4.5328e-08,-4.3900e-07,
                       9.7814e-12 ]

    return params

