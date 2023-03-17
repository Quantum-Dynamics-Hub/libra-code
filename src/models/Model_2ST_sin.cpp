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
  \file Model_2ST_sin.cpp
  \brief The file implements the 2-state sin-perturbed Rabi model potentials (different nuclear dimensions) 
*/

#include "Models_2_state.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


vector<double> set_params_2S_1D_sin(std::string model){
/**
  A helper function to setup the parameters used in the published works

*/

  vector<double> params(12, 0.0);   

  if(model=="SIN:Case1:JCC:2016:37:1626"){
    /**
    The sine-perturbed 2-state potential of Akimov
    [Akimov, A. V. Libra: An Open-Source “methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
     J. Comput. Chem. 2016, 37, 1626–1649]
    */                                                                               
    params[0] = -0.001;  params[1] = 0.0;    params[2] = 0.0;  params[3] = 1.0;  // H00
    params[4] = 0.001;   params[5] = 0.0;    params[6] = 0.0;  params[7] = 1.0;  // H01
    params[8] = 0.001;   params[9] = 0.0002; params[10]= 0.0;  params[11]= 1.0;  // H11

  }

  else if(model=="SIN:Case2:JCC:2016:37:1626"){
    /**
    The sine-perturbed 2-state potential of Akimov
    [Akimov, A. V. Libra: An Open-Source “methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
     J. Comput. Chem. 2016, 37, 1626–1649]
    */                                                                               
    params[0] = -0.001;  params[1] = 0.0;    params[2] = 0.0;  params[3] = 1.0;  // H00
    params[4] = 0.001;   params[5] = 0.0;    params[6] = 0.0;  params[7] = 1.0;  // H01
    params[8] = 0.001;   params[9] = 0.0019; params[10]= 0.0;  params[11]= 1.0;  // H11

  }

  else if(model=="SIN:Case3:JCC:2016:37:1626"){
    /**
    The sine-perturbed 2-state potential of Akimov
    [Akimov, A. V. Libra: An Open-Source “methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
     J. Comput. Chem. 2016, 37, 1626–1649]
    */                                                                               
    params[0] = -0.001;  params[1] = 0.0;    params[2] = 0.0;  params[3] = 1.0;  // H00
    params[4] = 0.001;   params[5] = 0.0;    params[6] = 0.0;  params[7] = 1.0;  // H01
    params[8] = 0.001;   params[9] = 0.0019; params[10]= 0.0;  params[11]= 5.0;  // H11

  }

  else if(model=="SIN:JPCL:2018:9:248"){
    /**
    Another sine-perturbed 2-state potential of Akimov.

    The parameters are chosen in units of thermal energy: E00 = -k_B * T, E11 = k_B * T, V = k_B * T,
    so that the minimal PES gap (splitting) is min [sqrt(H00-H11)^2 + 4H01^2)] = 2H01 = 2k_B *T, which 
    corresponds to the Boltzman factor of exp(-2) = 0.135. This magnitude ensures
    that the thermal excitations are minimized, but are still possible
    At the same time, min(H11) =0, max(H11) = 2k_B * T, and 1/2 [sqrt(H00-H11)^2 + 4H01^2)] > H01 = k_B * T
    So that the ground state PES barrier is smaller than the thermal energy, max(E0) - min(E0) < k_B * T
    Thus, the system can explore all the regions of the ground state PES under the influence of thermal fluctuations

    [Akimov, A. V. Stochastic and Quasi-Stochastic Hamiltonians for Long-Time Nonadiabatic Molecular Dynamics.
     J. Phys. Chem. Lett. 2017, 9, 248-257]

    */

    params[0] = -0.0009500434287;  params[1] = 0.0;             params[2] = 0.0;  params[3] = 1.0;  // H00
    params[4] = 0.0009500434287;   params[5] = 0.0;             params[6] = 0.0;  params[7] = 1.0;  // H01
    params[8] = 0.0009500434287;   params[9] = 0.0009500434287; params[10]= 0.0;  params[11]= 5.0;  // H11

  }

  return params;

}


void model_2S_1D_sin(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
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
   E0                      params[0]        -0.01
   A0                      params[1]         0.0
   x00                     params[2]         0.0
   L00                     params[3]         1.0

   V                       params[4]         0.005
   v                       params[5]         0.0
   x01                     params[6]         0.0
   L01                     params[7]         1.0

   E1                      params[8]         0.01
   A1                      params[9]         0.001
   x11                     params[10]        0.0
   L11                     params[11]        1.0

   
  Hamiltonian and its derivatives in diabatic representation:

       | E0 + A0 * sin(2*pi*(x-x00)/L00)     V + v * sin(2*pi*(x-x01)/L01)    |
   H = |                                                                      |
       |   V + v * sin(2*pi*(x-x01)/L01)     E1 + A1 * sin(2*pi*(x-x11)/L11)  |

   S = identity

   dc1 = 0.0

*/

  // Default parameters
  double E0 = -0.010;    double A0 = 0.0;    double x00 = 0.0;   double L00 = 1.0;  
  double V = 0.005;      double v = 0.0;     double x01 = 0.0;   double L01 = 1.0;  
  double E1 = 0.010;     double A1 = 0.001;  double x11 = 0.0;   double L11 = 1.0;  

  if(params.size()>=1){  E0 = params[0];    }
  if(params.size()>=2){  A0 = params[1];    }
  if(params.size()>=3){  x00= params[2];    }
  if(params.size()>=4){  L00= params[3];    }
  if(params.size()>=5){  V  = params[4];    }
  if(params.size()>=6){  v  = params[5];    }
  if(params.size()>=7){  x01= params[6];    }
  if(params.size()>=8){  L01= params[7];    }
  if(params.size()>=9){  E1 = params[8];    }
  if(params.size()>=10){ A1 = params[9];    }
  if(params.size()>=11){ x11= params[10];   }
  if(params.size()>=12){ L11= params[11];   }


  // H00 and dH00
  double argg = 2.0*M_PI/L00;
  double H00 = E0 + A0 * sin(argg*(q[0]-x00));
  double dH00 = A0 * argg * cos(argg *(q[0]-x00));

  // H01 and dH01
  argg = 2.0*M_PI/L01;
  double H01 = V + v * sin(argg*(q[0]-x01));
  double dH01 = v * argg * cos(argg *(q[0]-x01));

  // H11 and dH11
  argg = 2.0*M_PI/L11;
  double H11 = E1 + A1 * sin(argg*(q[0]-x11));
  double dH11 = A1 * argg * cos(argg *(q[0]-x11));


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



vector<double> set_params_2S_2D_sin(std::string model){
/**
  A helper function to setup the parameters used in the published works

*/

  vector<double> params(18, 0.0);   

  if(model=="SIN:JCC:2016:37:1626"){
    /**
    The sine-perturbed 2-state potential of Akimov
    [Akimov, A. V. Libra: An Open-Source “methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
     J. Comput. Chem. 2016, 37, 1626–1649]
    */                                                                               
    params[0]  =-0.001;  params[1]  = 0.0;    params[2]  = 0.0;  params[3]  = 1.0;  params[4]  = 0.0;  params[5]  = 1.0; // H00
    params[6]  = 0.001;  params[7]  = 0.0;    params[8]  = 0.0;  params[9]  = 1.0;  params[10] = 0.0;  params[11] = 1.0; // H01
    params[12] = 0.001;  params[13] = 0.0019; params[14] = 0.0;  params[15] = 1.0;  params[16] = 0.0;  params[17] = 1.0; // H11

  }

  return params;

}



void model_2S_2D_sin(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params){ 
/*** 
    To use with the nHamiltonian class

  \param[out] Hdia  The Hamiltonian in the diabatic basis (diabatic Hamiltonian)
  \param[out] Sdia  The overlap matrix in the diabatic basis
  \param[out] d1ham_dia  The 1-st order derivatives of the diabatic Hamiltonian w.r.t. all nuclear DOFs
  \param[out] dc1_dia  The 1-st order derivative couplings in the diabatic basis w.r.t. all nuclear DOFs
  \param[in] q The nuclear DOFs
  \param[in] params The model parameters: up to 18 parameters (see the chart below) will be used. If not defined,
            the default values will be used:

  Internal parameter        Input        Default value
   E0                      params[0]        -0.01
   A0                      params[1]         0.0
   x00                     params[2]         0.0
   L00_x                   params[3]         1.0
   y00                     params[4]         0.0
   L00_y                   params[5]         1.0

   V                       params[6]         0.005
   v                       params[7]         0.0
   x01                     params[8]         0.0
   L01_x                   params[9]         1.0
   y01                     params[10]        0.0
   L01_y                   params[11]        1.0

   E1                      params[12]        0.01
   A1                      params[13]        0.001
   x11                     params[14]        0.0
   L11_x                   params[15]        1.0
   y11                     params[16]        0.0
   L11_y                   params[17]        1.0


   
  Hamiltonian and its derivatives in diabatic representation:

       | E0 + A0 * sin(2*pi*(x-x00)/L00_x) * sin(2*pi*(y-y00)/L00_y)    V + v * sin(2*pi*(x-x01)/L01_x) * sin(2*pi*(y-y01)/L01_y)   |
   H = |                                                                                                                            |
       |   V + v * sin(2*pi*(x-x01)/L01_x) * sin(2*pi*(y-y01)/L01_y)   E1 + A1 * sin(2*pi*(x-x11)/L11_x) * sin(2*pi*(y-y11)/L11_y)  |

   S = identity

   dc1 = 0.0

*/

  // Default parameters
  double E0 = -0.010;    double A0 = 0.0;    double x00 = 0.0;   double L00_x = 1.0;  double y00 = 0.0;   double L00_y = 1.0;  
  double V = 0.005;      double v = 0.0;     double x01 = 0.0;   double L01_x = 1.0;  double y01 = 0.0;   double L01_y = 1.0;
  double E1 = 0.010;     double A1 = 0.001;  double x11 = 0.0;   double L11_x = 1.0;  double y11 = 0.0;   double L11_y = 1.0;

  if(params.size()>=1){  E0    = params[0];    }
  if(params.size()>=2){  A0    = params[1];    }
  if(params.size()>=3){  x00   = params[2];    }
  if(params.size()>=4){  L00_x = params[3];    }
  if(params.size()>=5){  y00   = params[4];    }
  if(params.size()>=6){  L00_y = params[5];    }
                            
  if(params.size()>=7){  V     = params[6];    }
  if(params.size()>=8){  v     = params[7];    }
  if(params.size()>=9){  x01   = params[8];    }
  if(params.size()>=10){ L01_x = params[9];    }
  if(params.size()>=11){ y01   = params[10];   }
  if(params.size()>=12){ L01_y = params[11];   }

  if(params.size()>=13){ E1    = params[12];   }
  if(params.size()>=14){ A1    = params[13];   }
  if(params.size()>=15){ x11   = params[14];   }
  if(params.size()>=16){ L11_x = params[15];   }
  if(params.size()>=17){ y11   = params[16];   }
  if(params.size()>=18){ L11_y = params[17];   }



  // H00 and dH00
  double omx = 2.0*M_PI/L00_x;    double sx = sin(omx*(q[0]-x00));
  double omy = 2.0*M_PI/L00_y;    double sy = sin(omy*(q[1]-y00));
  double H00 = E0 + A0 * sx * sy;
  double dH00x = A0 * omx * cos(omx *(q[0]-x00)) * sy;
  double dH00y = A0 * omy * sx * cos(omx *(q[1]-x00));


  // H01 and dH01
  omx = 2.0*M_PI/L01_x;    sx = sin(omx*(q[0]-x01));
  omy = 2.0*M_PI/L01_y;    sy = sin(omy*(q[1]-y01));
  double H01 = V + v * sx * sy;
  double dH01x = v * omx * cos(omx *(q[0]-x01)) * sy;
  double dH01y = v * omy * sx * cos(omx *(q[1]-x01));

  // H11 and dH11
  omx = 2.0*M_PI/L11_x;    sx = sin(omx*(q[0]-x11));
  omy = 2.0*M_PI/L11_y;    sy = sin(omy*(q[1]-y11));
  double H11 = E1 + A1 * sx * sy;
  double dH11x = A1 * omx * cos(omx *(q[0]-x11)) * sy;
  double dH11y = A1 * omy * sx * cos(omx *(q[1]-x11));



  Hdia.set(0,0, H00, 0.0);   Hdia.set(0,1, H01, 0.0); 
  Hdia.set(1,0, H01, 0.0);   Hdia.set(1,1, H11, 0.0); 

  Sdia.set(0,0, 1.0, 0.0);   Sdia.set(0,1, 0.0, 0.0); 
  Sdia.set(1,0, 0.0, 0.0);   Sdia.set(1,1, 1.0, 0.0); 

  //  d Hdia / dq_0
  d1ham_dia[0].set(0,0, dH00x, 0.0);  d1ham_dia[0].set(0,1, dH01x, 0.0);
  d1ham_dia[0].set(1,0, dH01x, 0.0);  d1ham_dia[0].set(1,1, dH11x, 0.0);

  //  d Hdia / dq_1
  d1ham_dia[1].set(0,0, dH00y, 0.0);  d1ham_dia[1].set(0,1, dH01y, 0.0);
  d1ham_dia[1].set(1,0, dH01y, 0.0);  d1ham_dia[1].set(1,1, dH11y, 0.0);


  //  <dia| d/dq_0| dia >
  dc1_dia[0].set(0,0, 0.0, 0.0);   dc1_dia[0].set(0,1, 0.0, 0.0); 
  dc1_dia[0].set(1,0, 0.0, 0.0);   dc1_dia[0].set(1,1, 0.0, 0.0); 

  //  <dia| d/dq_1| dia >
  dc1_dia[1].set(0,0, 0.0, 0.0);   dc1_dia[1].set(0,1, 0.0, 0.0); 
  dc1_dia[1].set(1,0, 0.0, 0.0);   dc1_dia[1].set(1,1, 0.0, 0.0); 


}




}// namespace libmodels
}// liblibra
