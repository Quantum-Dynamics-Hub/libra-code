/*********************************************************************************
* Copyright (C) 2018 Brendan A. Smith, Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_1ST_double_well.cpp
  \brief The file implements the 1-state double well potentials (different dimensions)   
*/

#include "Models_1_state.h"
#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_model namespace
namespace libhamiltonian_model{


void model_double_well_1S_1D(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                             vector<double>& q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 1 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   A                       params[0]         1.0

  double_well hamiltonian and its derivatives in diabatic representation:

   H_00 = A*(0.25x^4 - 0.5x^2)  
  dH_00 = A*(x - x^3)
*/

    // Default parameters
    double A = 1.0;

    if(params.size()>=1){   A = params[0];   }
    
    // H00
    double q2 = q[0]*q[0];
    double H00 = A*(0.25*q2*q2 - 0.5*q2); 
    double dH00 = A*q[0]*(q2 - 1.0);

    Hdia.set(0,0, H00, 0.0); 
    Sdia.set(0,0, 1.0, 0.0); 

    //  d Hdia / dq_0
    d1ham_dia[0].set(0,0, dH00, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0].set(0,0, 0.0, 0.0); 

}



}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra
