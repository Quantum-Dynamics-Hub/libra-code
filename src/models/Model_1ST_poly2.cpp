/*********************************************************************************
* Copyright (C) 2018-2022 Brendan A. Smith, Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_1ST_poly2.cpp
  \brief The file implements the 1-state harmonic potentials (different nuclear dimensions) = 
  = 2-nd order polynomials
*/

#include "Models_1_state.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


void model_1S_1D_poly2(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                       vector<double>& q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 2 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   k                       params[0]         1.0
   q0                      params[1]         0.0

  Hamiltonian and its derivatives in diabatic representation:

   H_00 = 0.5*k*(q-q0)^2
  dH_00 = k*(q - q0)
*/

    // Default parameters
    double k = 1.0;
    double q0 = 0.0;

    if(params.size()>=1){    k = params[0];    }
    if(params.size()>=2){   q0 = params[1];    }


    // H00
    double diff = q[0] - q0;
    double H00 = 0.5*k*diff * diff; 
    double dH00 = k*diff;

    Hdia.set(0,0, H00, 0.0); 
    Sdia.set(0,0, 1.0, 0.0); 

    //  d Hdia / dq_0
    d1ham_dia[0].set(0,0, dH00, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0].set(0,0, 0.0, 0.0); 

}


void model_1S_1D_poly2(CMATRIX* Hdia, CMATRIX* Sdia, vector<CMATRIX*>& d1ham_dia, vector<CMATRIX*>& dc1_dia,
                       vector<double>& q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 2 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   k                       params[0]         1.0
   q0                      params[1]         0.0

  Hamiltonian and its derivatives in diabatic representation:

   H_00 = 0.5*k*(q-q0)^2
  dH_00 = k*(q - q0)
*/

    // Default parameters
    double k = 1.0;
    double q0 = 0.0;

    if(params.size()>=1){    k = params[0];    }
    if(params.size()>=2){   q0 = params[1];    }


    // H00
    double diff = q[0] - q0;
    double H00 = 0.5*k*diff * diff; 
    double dH00 = k*diff;

    Hdia->set(0,0, H00, 0.0); 
    Sdia->set(0,0, 1.0, 0.0); 

    //  d Hdia / dq_0
    d1ham_dia[0]->set(0,0, dH00, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0]->set(0,0, 0.0, 0.0); 

}





}// namespace libmodels
}// liblibra

