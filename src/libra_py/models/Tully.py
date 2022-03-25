#*********************************************************************************                     
#* Copyright (C) 2018-2022 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: models_Tully
   :platform: Unix, Windows
   :synopsis: This module implements 3 Tully models for testing NA-MD dynamics
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

def Tully1_py(q, params, full_id):
    """
   
    Pure Python implementation of the Tully model I = Simple Avoided Crossing (SAC):

    H_00 = A*(1.0-exp(-B*x)) x>0,  
         = A*(exp(B*x)-1.0 ) x<0
    H_11 = -H_00
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.010, units: Ha]
            * **params["B"]** ( double ):  [ default: 1.600, units: Bohr^-1]
            * **params["C"]** ( double ):  [ default: 0.005, units: Ha]
            * **params["D"]** ( double ):  [ default: 1.000, units: Bohr^-2]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.010, "B":1.600, "C":0.005, "D":1.000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  
    #x = q.get(0)
    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    V11,dV11 = 0.0, 0.0
    if x>0:
        e = math.exp(-B*x)
        V11 = A*(1.0 - e)
        dV11 = A*B*e  
    else:
        e = math.exp(B*x)
        V11 = -A*(1.0 - e)
        dV11 = A*B*e

    e = math.exp(-D*x*x)
    V = C * e
    dV = -2.0*x*C*D*e

    Hdia.set(0,0, V11*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));      Hdia.set(1,1, V11*(-1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, dV11*(1.0+0.0j) );  d1ham_dia[i].set(0,1, dV*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dV*(1.0+0.0j));     d1ham_dia[i].set(1,1, dV11*(-1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def Tully1(q, params, full_id):
    """
   
    The implementation that calls the C++ implementation of Tully model I = Simple Avoided Crossing (SAC):

    H_00 = A*(1.0-exp(-B*x)) x>0,  
         = A*(exp(B*x)-1.0 ) x<0
    H_11 = -H_00
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.010, units: Ha]
            * **params["B"]** ( double ):  [ default: 1.600, units: Bohr^-1]
            * **params["C"]** ( double ):  [ default: 0.005, units: Ha]
            * **params["D"]** ( double ):  [ default: 1.000, units: Bohr^-2]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"A":0.010, "B":1.600, "C":0.005, "D":1.000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]


    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)

    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(x)
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C); prms.append(D);
    
    model_SAC(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj



def Tully2(q, params, full_id):
    """
   
    The implementation that calls the C++ implementation of Tully model II = Double Avoided Crossing (DAC):

    H_00 = 0.0
    H_11 = E - A*exp(-B*x^2)
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.100, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.028, units: Bohr^-2]
            * **params["C"]** ( double ):  [ default: 0.015, units: Ha]
            * **params["D"]** ( double ):  [ default: 0.060, units: Bohr^-2]
            * **params["E"]** ( double ):  [ default: 0.050, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"A":0.10, "B":0.28, "C":0.015, "D":0.060, "E":0.050 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    E = params["E"]



    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(x)
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C); prms.append(D);  prms.append(E);
    
    model_DAC(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj



def Tully3(q, params, full_id):
    """
   
    The implementation that calls the C++ implementation of Tully model III =
    Extended Coupling With Reflection  (ECWR):

    H_00 = A
    H_11 = -H_00
    H_01 = B*exp(C*x);          x <= 0
           B*(2.0 - exp(-C*x)); x > 0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.0006, units: Ha ]
            * **params["B"]** ( double ):  [ default: 0.1000, units: Ha ]
            * **params["C"]** ( double ):  [ default: 0.9000, units: Bohr^-1 ]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
 
    """

    critical_params = [ ] 
    default_params = {"A":0.0006, "B":0.1000, "C":0.9000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(x)
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C)
    
    model_ECWR(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj




def chain_potential(q, params, full_id):
    """    
    A 1D linear chain potential. 

    This is the Hamiltonian to mimic the one used by Tully and Parandekar to study 
    thermal equilibrium in quantum-classical systems

    Although the indended interpretation is that the first DOFs is quantum, this potential
    is totally generic, so all DOFs could be treated as quantum

    Args: 
        q ( MATRIX(ndof,1) ): coordinates of the nuclear DOFs
        params ( dictionary ): model parameters

            * **params["A"]** ( double ): [ units: None ]
            * **params["a"]** ( double ):  [ units: bohr-1 ]
            * **params["V0]** ( double ):  [ units: Ha ]            
            * **params["E_n"]** ( list of doubles ):  [ default: [0.0, 0.001, 0.001, 0.001], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 1.0, 1.0], units: Bohr]
            * **params["k_n"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001, 0.001], units: Ha/Bohr^2]
            * **params["V"]**   ( list of lists of double ):  [ default:  [[0.001]*4]*4, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(n,n) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(n,n) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(n, n) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(n, n) objects ): derivative coupling in the diabatic basis [ zero ]
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    X = q.col(indx)       # coordinates of all particles for this trajectory
    ndof = q.num_of_rows  # the total number of DOFs, the first one is the quantum DOF
    
    
    critical_params = ["E_n", "x_n", "k_n" ]
    default_params = { "V": [ [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001],
                              [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001] ],
                       "a":0.25, "U0":0.06665
                     }
    comn.check_input(params, default_params, critical_params)

    # Parameters describing the diabatic PES for quantum DOF
    E_n = params["E_n"] # energies of the minima of the diabatic PESs for quantum particle
    x_n = params["x_n"] # positions of diabatic PESs
    k_n = params["k_n"] # curvatures of diabatic PESs
    V = params["V"]     # diabatic couplings
    nstates = len(E_n)  # the number of electronic states
    
    # Parameters describing classical particles and their coupling to the quantum particle    
    a  = params["a"]
    U0 = params["U0"]
    a2 =  U0*a**2
    a3 = -U0*a**3
    a4 = 0.58*U0*a**4
    
                
    Hdia = CMATRIX(nstates,nstates)
    Sdia = CMATRIX(nstates,nstates)
    Sdia.identity()
    basis_transform = CMATRIX(nstates,nstates)
    basis_transform.identity()
    
    
    d1ham_dia = CMATRIXList();
    dc1_dia   = CMATRIXList();
    for k in range(ndof):
        d1ham_dia.append( CMATRIX(nstates,nstates) )
        dc1_dia.append( CMATRIX(nstates,nstates) ) 

        
    #========= Classical (bath) contributions =========    
    H_bath  = 0.0
    dH_bath = []
    for k in range(ndof-1):
        x1 = X.get(k)
        x2 = X.get(k+1)
        d = x1 - x2
        d2 = d*d
        d3 = d2*d
        d4 = d3*d
        H_bath += (a2*d2 + a3*d3 + a4*d4)
        dH_bath.append( 2.0*a2*d + 3.0*a3*d2 + 4.0*a4*d3 )

        
    for i in range(nstates):        
        Hdia.add(i,i, H_bath*(1.0+0.0j))
                
        k = 0
        d1ham_dia[k].add(i, i, ( dH_bath[k] ) * (1.0+0.0j) )
        
        for k in range(1, ndof-1):
            d1ham_dia[k].add(i, i, ( dH_bath[k] - dH_bath[k-1] ) * (1.0+0.0j) )
            
        k = ndof - 1
        d1ham_dia[k].add(i, i, ( -dH_bath[k-1] ) * (1.0+0.0j) )
        
    
        
    """
    for k in range(1, ndof):        
        
        d = (X.get(0) - X.get(k) )        
        d2 = d*d        
        H_bath = 0.5 * U0*d2
        dH_bath = U0 * d
        
        
        for i in range(nstates):        
            Hdia.add(i,i, H_bath*(1.0+0.0j) )            
            d1ham_dia[k].add(i, i, ( -dH_bath ) * (1.0+0.0j) )
            d1ham_dia[0].add(i, i, (  dH_bath ) * (1.0+0.0j) )
        
    """
        
    #========== Quantum system ==================    
    x = X.get(0)    
    for i in range(nstates):        
        # Energies
        Hdia.add(i,i,  (E_n[i] + 0.5*k_n[i]*(x - x_n[i])**2) * (1.0+0.0j) )
                
        for j in range(nstates):
            if i!=j:
                Hdia.add(i,j,  V[i][j] * (1.0+0.0j) )
        
        # Derivatives
        k = 0 # quantum DOF
        d1ham_dia[k].add(i,i, (k_n[i] * (x - x_n[i]))*(1.0+0.0j) )

                        
    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia
    #obj.basis_transform = basis_transform
    
    return obj

