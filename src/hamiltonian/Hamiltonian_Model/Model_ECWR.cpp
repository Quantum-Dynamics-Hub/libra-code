/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Model_ECWR.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libhamiltonian{
namespace libhamiltonian_model{


void ECWR_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params){ 
// ECWR hamiltonian in diabatic representation

  if(H->n_elts!=4){ std::cout<<"Error in ECWR_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH->n_elts!=4){ std::cout<<"Error in ECWR_Ham: dH matrix must be allocated\n"; exit(0);}
  if(d2H->n_elts!=4){ std::cout<<"Error in ECWR_Ham: d2H matrix must be allocated\n"; exit(0);}

  // Default parameters
  double A = 0.0006;  double B = 0.100;
  double C = 0.900;   double D = 1.00;  

  if(params.size()>=4){
    A = params[0];    B = params[1];
    C = params[2];    D = params[3];
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

}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra
