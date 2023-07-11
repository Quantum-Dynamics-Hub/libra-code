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

#include "Model_SEXCH.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libmodels namespace
namespace libmodels{


void SEXCH_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
// Superexchange Hamiltonian in diabatic representation

  if(H->n_elts!=9){ std::cout<<"Error in SEXCH_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->n_elts!=9){ std::cout<<"Error in SEXCH_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->n_elts!=9){ std::cout<<"Error in SEXCH_Ham: d2H matrix must be allocated\n"; exit(0);}


  // Default parameters
  double A1=0.001;  double A2=0.010;  double A3=0.000;
  double B1=0.500;  double B2=0.500;  double B3=0.500;
  double C1=0.000;  double C2=0.000;  double C3=0.000;
  double D1=0.000;  double D2=0.010;  double D3=0.005;

  if(params.size()>=12){
    A1 = params[0];   A2 = params[1];   A3 = params[2];
    B1 = params[3];   B2 = params[4];   B3 = params[5];
    C1 = params[6];   C2 = params[7];   C3 = params[8];
    D1 = params[9];   D2 = params[10];  D3 = params[11];
  }


  // H 
  H->M[0] = D1;          H->M[1] = A1*exp(-B1*(x-C1)*(x-C1));   H->M[2] = A3*exp(-B3*(x-C3)*(x-C3)); 
  H->M[3] = H->M[1];     H->M[4] = D2;                          H->M[5] = A2*exp(-B2*(x-C2)*(x-C2));  
  H->M[6] = H->M[2];     H->M[7] = H->M[5];                     H->M[8] = D3;

  // dH/dx
  dH->M[0] = 0.0;        dH->M[1] = -2.0*B1*(x-C1)*H->M[1];     dH->M[2] = -2.0*B3*(x-C3)*H->M[2]; 
  dH->M[3] = dH->M[1];   dH->M[4] = 0.0;                        dH->M[5] = -2.0*B2*(x-C2)*H->M[5];
  dH->M[6] = dH->M[2];   dH->M[7] = dH->M[5];                   dH->M[8] = 0.0;

  // d2H/dx2
  d2H->M[0] = 0.0;       d2H->M[1] = -2.0*B1*(H->M[1] + (x-C1)*dH->M[1]);  d2H->M[2] = -2.0*B3*(H->M[2] + (x-C3)*dH->M[2]); 
  d2H->M[3] = dH->M[1];  d2H->M[4] = 0.0;                                  d2H->M[5] = -2.0*B2*(H->M[5] + (x-C2)*dH->M[5]);
  d2H->M[6] = dH->M[2];  d2H->M[7] = dH->M[5];                             d2H->M[8] = 0.0;


}

boost::python::list SEXCH_Ham(double x, boost::python::list params_){ 
// SEXCH hamiltonian in diabatic representation

  MATRIX H(3,3);
  MATRIX dH(3,3);
  MATRIX d2H(3,3);

  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  SEXCH_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libmodels
}// liblibra

