/*********************************************************************************
* Copyright (C) 2017-2022 Brendan A. Smith , Alexey V. Akimov
* Copyright (C) 2015-2017 Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_double_well.cpp
  \brief The file implements the functions for computing the dobule_well Hamiltonian and its derivatives - model for 
  tunneling between two symmetric potentials seperated by a barrier.
    
*/

#include "Model_double_well.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


void double_well_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
/** The symmetric double well potential. Two local minima are located at x = -1.0 and 1.0. A local maximum is located at x = 0.

  \param[in] x The nuclear coordinate (1D)
  \param[out] H The pointer to the matrix in which the Hamiltonian will be written
  \param[out] dH The pointer to the matrix in which the 1-st order derivatives of the Hamiltonian will be written
  \param[out] d2H The pointer to the matrix in which the 2-nd order derivatives of the Hamiltonian will be written
  \param[in] Currently, 0 parameters are used for the double well potential.
  values will be used:

  Double well potential Hamiltonian and its derivatives in diabatic representation:
  
  H_00 = A * ( 0.25*x^4 - 0.5*x^2 )
 
  ### explicit functions ###

  f(x) = A * [ (1/4)*x^4 - (1/2)*x^2 ]
  f'(x) = A * [ x^3 - x ]
  f''(x) = A * [ 3x^2 - 1 ]
*/


  if(H->n_elts!=1){ std::cout<<"Error in double_well_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->n_elts!=1){ std::cout<<"Error in double_well_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->n_elts!=1){ std::cout<<"Error in double_well_Ham: d2H matrix must be allocated\n"; exit(0);}


  // Default parameters: 
  double A = 1.0; //< sets depth of wells to 0.25 atomic energy units. 
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;

  if(params.size()>=1){   A = params[0];   }

  // H00
  H->M[0] = A * (0.25*x4 - 0.5*x2);
  dH->M[0] = A * (x3 - x);
  d2H->M[0] = A * (3*x2 - 1);
  
}

boost::python::list double_well_Ham(double x, boost::python::list params_){ 
/** The symmetric double_well potential 

  This is Python-friendly version

  \param[in] x The nuclear coordinate (1D)
  \param[out] H The pointer to the matrix in which the Hamiltonian will be written
  \param[out] dH The pointer to the matrix in which the 1-st order derivatives of the Hamiltonian will be written
  \param[out] d2H The pointer to the matrix in which the 2-nd order derivatives of the Hamiltonian will be written
  \param[in] Currently, 0 parameters are used for the double well potential.
  values will be used:

  double_well hamiltonian and its derivatives in diabatic representation:

  H_00 = A * [ 0.25*x^4 - 0.5*x^2 ]

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

  double_well_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libmodels
}// liblibra
