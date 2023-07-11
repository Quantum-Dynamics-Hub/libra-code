/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_cubic.cpp
  \brief The file implements the functions for computing cubic Hamiltonian and its derivatives - model for 
  decay from a metastable state
  see Donoso, A.; Martens, C. Quantum Tunneling Using Entangled Classical Trajectories. Phys. Rev. Lett. 2001, 87, 223202.
    
*/

#include "Model_cubic.h"

/// liblibra namespace
namespace liblibra{

/// libmodels namespace
namespace libmodels{


void cubic_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
/** The metastable potential - e.g. see Donoso, A.; Martens, C. Quantum Tunneling Using Entangled Classical Trajectories. Phys. Rev. Lett. 2001, 87, 223202.
  The default parameters are chosen to make 2 bound state, so the system is highly quantum-mechanical

  \param[in] x The nuclear coordinate (1D)
  \param[out] H The pointer to the matrix in which the Hamiltonian will be written
  \param[out] dH The pointer to the matrix in which the 1-st order derivatives of the Hamiltonian will be written
  \param[out] d2H The pointer to the matrix in which the 2-nd order derivatives of the Hamiltonian will be written
  \param[in] params The model parameters: can be up to 3 parameters (see the chart below). Otherwise, the default
  values will be used:

  Internal parameter        Input        Default value
   m                       param[0]         2000.0
   w                       param[1]         0.010
   b                       param[2]         0.2981

  Cubic hamiltonian and its derivatives in diabatic representation:

  H_00 = 0.5*m*w^2 * x^2 - (b/3) * x^3
*/


  if(H->n_elts!=1){ std::cout<<"Error in cubic_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->n_elts!=1){ std::cout<<"Error in cubic_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->n_elts!=1){ std::cout<<"Error in cubic_Ham: d2H matrix must be allocated\n"; exit(0);}

  double e;

  // Default parameters: this corresponds to the barrier height: V_Eq = 0.015   at   q_eq = 0.6709
  double m = 2000.0;  ///< mass m
  double w = 0.010;   ///< omega 
  double b = 0.2981;  ///< cubic term
  double mw2 = m*w*w;
  double x2 = x*x;
  double x3 = x2*x;
  

  if(params.size()>=1){   m = params[0];    }
  if(params.size()>=2){   w = params[1];    }
  if(params.size()>=3){   b = params[2];    }


  // H00
  H->M[0] = 0.5*mw2*x2 - (b/3.0)*x3;
  dH->M[0] = mw2*x - b*x2;
  d2H->M[0] = mw2 - 2.0*b*x;

  
}

boost::python::list cubic_Ham(double x, boost::python::list params_){ 
/** The metastable potential - e.g. see Donoso, A.; Martens, C. Quantum Tunneling Using Entangled Classical Trajectories. Phys. Rev. Lett. 2001, 87, 223202.
  The default parameters are chosen to make 2 bound state, so the system is highly quantum-mechanical

  This is Python-friendly version

  \param[in] x The nuclear coordinate (1D)
  \param[out] H The pointer to the matrix in which the Hamiltonian will be written
  \param[out] dH The pointer to the matrix in which the 1-st order derivatives of the Hamiltonian will be written
  \param[out] d2H The pointer to the matrix in which the 2-nd order derivatives of the Hamiltonian will be written
  \param[in] params The model parameters: can be up to 3 parameters (see the chart below). Otherwise, the default
  values will be used:

  Internal parameter        Input        Default value
   m                       param[0]         2000.0
   w                       param[1]         0.010
   b                       param[2]         0.2981

  Cubic hamiltonian and its derivatives in diabatic representation:

  H_00 = 0.5*m*w^2 * x^2 - (b/3) * x^3

  Returns the Python list, res, of the following objects:
  res[0] = x - coordinate, 
  res[1] = H - the Hamiltonian matrix 
  res[2] = dH - the Hamiltonian derivatives of the 1-st order
  res[3] = d2H - the Hamiltonian derivatives of the 2-nd order

*/

  MATRIX H(2,2);
  MATRIX dH(2,2);
  MATRIX d2H(2,2);

  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  cubic_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libmodels
}// liblibra
