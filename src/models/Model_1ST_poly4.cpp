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
  \file Model_1ST_poly4.cpp
  \brief The file implements the 1-state anharmonic potentials (different nuclear dimensions) = 
  = 4-th order-polynomials
*/

#include "Models_1_state.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


vector<double> set_params_1S_1D_poly4(std::string model){
/**
  A helper function to setup the parameters used in the published works

*/

  vector<double> params(4, 0.0);   

  if(model=="Cubic:PRL:2001:87:223202"){
    /**
    The metastable state potential of Donoso and Martens
    [Donoso, A.; Martens, C. Quantum Tunneling Using Entangled Classical Trajectories. Phys. Rev. Lett. 2001, 87, 223202]
    These parameters will yield 2 bound state, so the system is highly quantum-mechanical
    */
  
    params[0] = 0.2; 
    params[1] = 0.0;
    params[2] = -0.09936667;
    params[3] = 0.0;
  }

  else if(model=="Cubic:JCP:2000:113:6557"){
    /**
    The metastable state potential of Prezhdo and Pereverzev
    [Prezhdo, O. V.; Pereverzev, Y. V. Quantized Hamilton Dynamics. J. Chem. Phys. 2000, 113, 6557–6565]
    */
    params[0] = 1.0; 
    params[1] = 0.0;
    params[2] = 0.1;
    params[3] = 0.0;

  }

  else if(model=="DoubleWell:JCP:2002:117:2995"){
    /**
    The double well potential of Prezhdo
    [Prezhdo, O. V. Classical Mapping for Second-Order Quantized Hamiltonian Dynamics. J. Chem. Phys. 2002, 117, 2995–3002]
    */
    params[0] =-5.0; 
    params[1] = 0.0;
    params[2] = 0.0;
    params[3] = 1.0;

  }

  else if(model=="AHO:"){
    /**
    The anharmonic oscillator of Ananth [??]
    */

    params[0] = 1.0; 
    params[1] = 0.0;
    params[2] = -0.1;
    params[3] = 0.1;

  }

  else if(model=="DW:"){ 
    /**
    For the double well potential of Smith [under review]
    */

    params[0] = -0.5; 
    params[1] = 0.0;
    params[2] = 0.0;
    params[3] = 0.25;

  }

  return params;

}


void model_1S_1D_poly4(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                       vector<double>& q, vector<double>& params){ 
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
   k                       params[0]         1.0
   q0                      params[1]         0.0
   a                       params[2]         0.0
   b                       params[3]         0.0
   
  Hamiltonian and its derivatives in diabatic representation:

   H_00 = 0.5*k*(q-q0)^2 + a * (q-q0)^3 + b * (q-q0)^4
  dH_00 = k*(q - q0) + 3*a * (q-q0)^2 + 4*b * (q-q0)^3

*/


    // Default parameters
    double k = 1.0; // harmonic constant
    double q0= 0.0; // center of the potential
    double a = 0.0; // cubic correction
    double b = 0.0; // quartic correction

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

void model_1S_1D_poly4(CMATRIX* Hdia, CMATRIX* Sdia, vector<CMATRIX*>& d1ham_dia, vector<CMATRIX*>& dc1_dia,
                       vector<double>& q, vector<double>& params){ 
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
   k                       params[0]         1.0
   q0                      params[1]         0.0
   a                       params[2]         0.0
   b                       params[3]         0.0
   
  Hamiltonian and its derivatives in diabatic representation:

   H_00 = 0.5*k*(q-q0)^2 + a * (q-q0)^3 + b * (q-q0)^4
  dH_00 = k*(q - q0) + 3*a * (q-q0)^2 + 4*b * (q-q0)^3

*/


    // Default parameters
    double k = 1.0; // harmonic constant
    double q0= 0.0; // center of the potential
    double a = 0.0; // cubic correction
    double b = 0.0; // quartic correction

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

    Hdia->set(0,0, H00, 0.0); 
    Sdia->set(0,0, 1.0, 0.0); 

    //  d Hdia / dq_0
    d1ham_dia[0]->set(0,0, dH00, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0]->set(0,0, 0.0, 0.0); 

}





}// namespace libmodels
}// liblibra
