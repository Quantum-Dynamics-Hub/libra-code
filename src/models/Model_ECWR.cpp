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

#include "Model_ECWR.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libmodels namespace
namespace libmodels{


void model_ECWR(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
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
   A                       params[0]         0.0006
   B                       params[1]         0.1000
   C                       params[2]         0.9000

  ECWR hamiltonian and its derivatives in diabatic representation:

  H_00 = A
  H_11 = -H_00
  H_01 = B*exp(C*x);          x <= 0
         B*(2.0 - exp(-C*x)); x > 0
*/

    double e;

    // ECWR potetnial
    // Default parameters
    double A = 0.0006;  double B = 0.100;   double C = 0.900;  

    if(params.size()>=3){
      A = params[0];    B = params[1];     C = params[2];   
    }

    double H01, dH01;

    // H01
    if(q[0]>=0){  e = exp(-C*q[0]);  H01 = B*(2.0 - e);   dH01 = B*C*e;  } 
    else{      e = exp(C*q[0]);  H01 = B*e;  dH01 = B*C*e;   }


    Sdia.set(0,0, 1.0, 0.0);  Sdia.set(0,1, 0.0, 0.0);
    Sdia.set(1,0, 0.0, 0.0);  Sdia.set(1,1, 1.0, 0.0);

    Hdia.set(0,0, A, 0.0);    Hdia.set(0,1,  H01, 0.0);
    Hdia.set(1,0, H01, 0.0);  Hdia.set(1,1, -A, 0.0);

    //  d Hdia / dq_0
    d1ham_dia[0].set(0,0, 0.0, 0.0);   d1ham_dia[0].set(0,1, dH01, 0.0);
    d1ham_dia[0].set(1,0, dH01, 0.0);   d1ham_dia[0].set(1,1, 0.0, 0.0);

    //  <dia| d/dq_0| dia >
    dc1_dia[0].set(0,0, 0.0, 0.0);   dc1_dia[0].set(0,1, 0.0, 0.0);
    dc1_dia[0].set(1,0, 0.0, 0.0);   dc1_dia[0].set(1,1, 0.0, 0.0);


}




void ECWR_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
// ECWR hamiltonian in diabatic representation

  if(H->n_elts!=4){ std::cout<<"Error in ECWR_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->n_elts!=4){ std::cout<<"Error in ECWR_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->n_elts!=4){ std::cout<<"Error in ECWR_Ham: d2H matrix must be allocated\n"; exit(0);}

  // Default parameters
  double A = 0.0006;  double B = 0.100;  double C = 0.900;   

  if(params.size()>=3){
    A = params[0];    B = params[1];    C = params[2];   
  }


  // H00
  H->M[0] =  A;    
  dH->M[0] = 0.0;
  d2H->M[0] = 0.0;

  // H11
  H->M[3] = -A;  
  dH->M[3] = 0.0;
  d2H->M[3] = 0.0;


  // H01 and H10
  if(x<=0){ 
    H->M[1] = H->M[2] = B*exp(C*x);
    dH->M[1] = dH->M[2] = C*H->M[1]; 
    d2H->M[1] = d2H->M[2] = C*dH->M[1]; 
  }
  else{ 
    double e = exp(-C*x); 
    H->M[1] = H->M[2] = B*(2.0 - e);
    dH->M[1] = dH->M[2] = B*C*e;
    d2H->M[1] = d2H->M[2] = -C*dH->M[1];
  }


}

boost::python::list ECWR_Ham(double x, boost::python::list params_){ 
// ECWR hamiltonian in diabatic representation

  MATRIX H(2,2);
  MATRIX dH(2,2);
  MATRIX d2H(2,2);

  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  ECWR_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libmodels
}// liblibra
