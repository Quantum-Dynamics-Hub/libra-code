/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Model_SAC.h"

namespace libhamiltonian{
namespace libhamiltonian_model{


void SAC_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
// SAC hamiltonian in diabatic representation

  if(H->num_of_elems!=4){ std::cout<<"Error in SAC_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->num_of_elems!=4){ std::cout<<"Error in SAC_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->num_of_elems!=4){ std::cout<<"Error in SAC_Ham: d2H matrix must be allocated\n"; exit(0);}


  double e;

  // SAC potetnial
  // Default parameters
  double A = 0.010;  double B = 1.600;
  double C = 0.005;  double D = 1.00;

  if(params.size()>=4){
    A = params[0];    B = params[1];
    C = params[2];    D = params[3];
  }

  // H00
  if(x>=0){ 
    e = exp(-B*x); 
    H->M[0] = A*(1.0 - e);
    dH->M[0] = A*B*e;
    d2H->M[0] = -B*dH->M[0];
  } 
  else{
    e = exp(B*x);
    H->M[0] = A*(e - 1.0);
    dH->M[0] = A*B*e;
    d2H->M[0] = B*dH->M[0];
  }

  // H11
  H->M[3] = -H->M[0]; 
  dH->M[3] = -dH->M[0];
  d2H->M[3] = -d2H->M[0];
  
  // H01 and H10
  H->M[1] = H->M[2] = C*exp(-D*x*x);
  dH->M[1] = dH->M[2] = -2.0*D*x*H->M[1];
  d2H->M[1] = d2H->M[2] = (-2.0*D*H->M[1] - 2.0*D*x*dH->M[1]);


}

boost::python::list SAC_Ham(double x, boost::python::list params_){ 
// SAC hamiltonian in diabatic representation

  MATRIX H(2,2);
  MATRIX dH(2,2);
  MATRIX d2H(2,2);

  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  SAC_Ham(x,&H,&dH,&d2H,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH);
  res.append(d2H);
 
  return res;

}

}// namespace libhamiltonian_model
}// namespace libhamiltonian

