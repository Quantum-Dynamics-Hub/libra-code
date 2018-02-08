/*********************************************************************************
* Copyright (C) 2015-2017 Brendan A. Smith , Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
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
#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_model namespace
namespace libhamiltonian_model{

void model_double_well(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia, vector<double> q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 4 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   A                       params[0]         1.0

  double_well hamiltonian and its derivatives in diabatic representation:

   H_00 = A*(0.25x^4 - 0.5x^2)  
  dH_00 = A*(x - x^3)
*/

    // symmetric_double_well potetnial
    // Default parameters
    double A = 1.0;

    if(params.size()>=1){
      A = params[0];
    }

    double H00, H01;
    double dH00, dH01;

    // H00;  H11 = 0.0 (single PES)
    H00 = A*(0.25*q[0]*q[0]*q[0]*q[0] - 0.5*q[0]*q[0]);  dH00 = A*(q[0]*q[0]*q[0] - q[0]);
    double H11 = 0.0; double dH11 = 0.0;

    // H01 = H10
    H01 = 0.0;   dH01 = 0.0;

    Sdia.set(0,0, 1.0, 0.0);  Sdia.set(0,1, 0.0, 0.0);
    Sdia.set(1,0, 0.0, 0.0);  Sdia.set(1,1, 1.0, 0.0);

    Hdia.set(0,0, H00, 0.0);  Hdia.set(0,1,  H01, 0.0);
    Hdia.set(1,0, H01, 0.0);  Hdia.set(1,1,  H11, 0.0);

    //  d Hdia / dq_0
    d1ham_dia[0].set(0,0, dH00, 0.0);   d1ham_dia[0].set(0,1, dH01, 0.0);
    d1ham_dia[0].set(1,0, dH01, 0.0);   d1ham_dia[0].set(1,1, dH11, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0].set(0,0, 0.0, 0.0);   dc1_dia[0].set(0,1, 0.0, 0.0);
    dc1_dia[0].set(1,0, 0.0, 0.0);   dc1_dia[0].set(1,1, 0.0, 0.0);

}

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

}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra
