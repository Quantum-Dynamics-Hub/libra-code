/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_2ST_tanh.cpp
  \brief The file implements the 2-state tanh-barrier potentials (different nuclear dimensions) 
*/

#include "Models_2_state.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


vector<double> set_params_2S_1D_tanh(std::string model){
/**
  A helper function to setup the parameters used in the published works

*/

  vector<double> params(7, 0.0);   

  if(model=="TANH:Case1:JCP:2018:148:102326"){
    /**
    The tangent barrier of Church (Ananth group)
    [Church, M. S.; Hele, T. J. H.; Ezra, G. S.; Ananth, N. Nonadiabatic Semiclassical Dynamics
     in the Mixed Quantum-Classical Initial Value Representation. The Journal of Chemical Physics 2018, 148, 102326]

    Similar setup is here:
    [Ananth, N.; Venkataraman, C.; Miller, W. H. Semiclassical Description of Electronically
     Nonadiabatic Dynamics via the Initial Value Representation. The Journal of Chemical Physics 2007, 127, 084114.]
    */                                                                               
    params[0] = 0.01;  params[1] = 1.6;  // H00
    params[2] = 0.01;  params[3] = 1.6;  // H11
    params[4] = 0.005; params[5] = 1.0;  params[6] = 0.0;  // H01

  }

  else if(model=="TANH:Case2:JCP:2018:148:102326"){
    /**
    The tangent barrier of Church (Ananth group)
    [Church, M. S.; Hele, T. J. H.; Ezra, G. S.; Ananth, N. Nonadiabatic Semiclassical Dynamics
     in the Mixed Quantum-Classical Initial Value Representation. The Journal of Chemical Physics 2018, 148, 102326]

    Similar setup is here:
    [Ananth, N.; Venkataraman, C.; Miller, W. H. Semiclassical Description of Electronically
     Nonadiabatic Dynamics via the Initial Value Representation. The Journal of Chemical Physics 2007, 127, 084114.]
    */                                                                               
    params[0] = 0.04;  params[1] = 1.0;  // H00
    params[2] = 0.01;  params[3] = 1.0;  // H11
    params[4] = 0.005; params[5] = 1.0;  params[6] = 0.7;  // H01

  }

  return params;

}


void model_2S_1D_tanh(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                      vector<double>& q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 7 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   V0                      params[0]         0.01
   alp0                    params[1]         1.0
   V1                      params[2]         0.01
   alp1                    params[3]         1.0
   a                       params[4]         0.005
   b                       params[5]         1.0
   x0                      params[6]         0.0

   
  Hamiltonian and its derivatives in diabatic representation:

       |   V0 * (1 + tanh(alp1 * x) )     a * exp(-b*(x + x0)^2 )    |
   H = |                                                             |
       |   a * exp(-b*(x + x0)^2 )      V1 * (1 - tanh(alp2 * x) )   |

   S = identity

   dc1 = 0.0

*/

  // Default parameters
  double V0 = 0.01;    double alp0 = 1.0;  
  double V1 = 0.01;    double alp1 = 1.0;  
  double a = 0.005;    double b = 1.0;    double x0 = 0.0;

  if(params.size()>=1){  V0   = params[0];    }
  if(params.size()>=2){  alp0 = params[1];    }
  if(params.size()>=3){  V1   = params[2];    }
  if(params.size()>=4){  alp1 = params[3];    }
  if(params.size()>=5){  a    = params[4];    }
  if(params.size()>=6){  b    = params[5];    }
  if(params.size()>=7){  x0   = params[6];    }


  // H00 and dH00
  double ta = tanh(alp0*q[0]);
  double H00 = V0 * (1.0 + ta);
  double dH00 = V0 * alp0 * (1.0 - ta*ta);

  // H11 and dH11
  ta = tanh(alp1*q[0]);
  double H11 = V1 * (1.0 - ta);
  double dH11 = -V1 * alp1 * (1.0 - ta*ta);

  // H01 and dH01
  double e = exp(-b*(q[0]+x0)*(q[0]+x0));
  double H01 = a * e;
  double dH01 = -2.0 * a * b * (q[0]+x0) * e;



  Hdia.set(0,0, H00, 0.0);   Hdia.set(0,1, H01, 0.0); 
  Hdia.set(1,0, H01, 0.0);   Hdia.set(1,1, H11, 0.0); 

  Sdia.set(0,0, 1.0, 0.0);   Sdia.set(0,1, 0.0, 0.0); 
  Sdia.set(1,0, 0.0, 0.0);   Sdia.set(1,1, 1.0, 0.0); 

  //  d Hdia / dq_0
  d1ham_dia[0].set(0,0, dH00, 0.0);  d1ham_dia[0].set(0,1, dH01, 0.0);
  d1ham_dia[0].set(1,0, dH01, 0.0);  d1ham_dia[0].set(1,1, dH11, 0.0);

  //  <dia| d/dq_0| dia >
  dc1_dia[0].set(0,0, 0.0, 0.0);   dc1_dia[0].set(0,1, 0.0, 0.0); 
  dc1_dia[0].set(1,0, 0.0, 0.0);   dc1_dia[0].set(1,1, 0.0, 0.0); 


}



}// namespace libmodels
}// liblibra
