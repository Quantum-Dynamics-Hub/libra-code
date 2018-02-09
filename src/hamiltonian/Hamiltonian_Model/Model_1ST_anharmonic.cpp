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
  \file Model_1ST_anharmonic.cpp
  \brief The file implements the 1-state anharmonic potentials (different nuclear dimensions)
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


void model_anharmonic_1S_1D(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
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
   a                       params[2]        -0.1
   b                       params[3]         0.1

  double_well hamiltonian and its derivatives in diabatic representation:

   H_00 = 0.5*k*(q-q0)^2 + a * (q-q0)^3 + b * (q-q0)^4
  dH_00 = k*(q - q0) + 3*a * (q-q0)^2 + 4*b * (q-q0)^3
*/


    // Default parameters
    double k = 1.0;  // harmonic constant
    double q0 = 0.0; // center of the potential
    double a = -0.1; // cubic correction
    double b = 0.1;  // quartic correction

    if(params.size()>=1){    k = params[0];    }
    if(params.size()>=2){   q0 = params[1];    }
    if(params.size()>=3){   a  = params[2];    }
    if(params.size()>=4){   b  = params[3];    }


    // H00
    double d = q[0] - q0;
    double d2 = d*d;
    double d3 = d2 * d;

    double H00 = 0.5*k*d2 + a*d3 + b*d3*d; 
    double dH00 = k*d - 3*a*d2 + 4*b*d3;

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
