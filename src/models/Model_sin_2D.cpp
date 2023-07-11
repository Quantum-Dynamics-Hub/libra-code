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

#include "Model_sin.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


void sin_2D_Ham(double x, double y, MATRIX* H, 
                MATRIX* dH1,  MATRIX* dH2, 
                MATRIX* d2H1, MATRIX* d2H2, vector<double>& params){ 
// sin_2D hamiltonian in diabatic representation
// dH1 - first derivative of H w.r.t. x
// dH2 - first derivative of H w.r.t. y
// d2H1 - second derivative of H w.r.t. x
// d2H2 - second derivative of H w.r.t. y
// no cross-derivatives yet


  if(H->n_elts!=4){ std::cout<<"Error in sin_2D_Ham: H matrix must be allocated\n"; exit(0);}
  if(dH1->n_elts!=4){ std::cout<<"Error in sin_2D_Ham: dH1 matrix must be allocated\n"; exit(0);}
  if(d2H1->n_elts!=4){ std::cout<<"Error in sin_2D_Ham: d2H1 matrix must be allocated\n"; exit(0);}
  if(dH2->n_elts!=4){ std::cout<<"Error in sin_2D_Ham: dH2 matrix must be allocated\n"; exit(0);}
  if(d2H2->n_elts!=4){ std::cout<<"Error in sin_2D_Ham: d2H2 matrix must be allocated\n"; exit(0);}


  double e;


  // Default parameters
  double A = -0.010;  double B = 0.010;
  double C = 0.005;   double D = 0.001;
  double Lx = 1.0;    double Ly = 1.0;


  if(params.size()>=6){
    A = params[0];    B = params[1];
    C = params[2];    D = params[3];
    Lx = params[4];   Ly = params[5];
  }

  // H00
  H->M[0] = A;
  dH1->M[0] = 0.0;
  d2H1->M[0] = 0.0;
  dH2->M[0] = 0.0;
  d2H2->M[0] = 0.0;

  // H11
  H->M[3] = B; 
  dH1->M[3] = 0.0;
  d2H1->M[3] = 0.0;
  dH2->M[3] = 0.0;
  d2H2->M[3] = 0.0;


  // H01 and H10
  H->M[1] = H->M[2] = C;
  dH1->M[1] = dH1->M[2] = 0.0;
  d2H1->M[1] = d2H1->M[2] = 0.0;
  dH2->M[1] = dH2->M[2] = 0.0;
  d2H2->M[1] = d2H2->M[2] = 0.0;


  
  // Now the correction to E1
  double om_x = 2.0*(M_PI/Lx);
  double om_y = 2.0*(M_PI/Ly);

  double cs_x = cos(om_x*x);
  double cs_y = cos(om_y*y);
  double si_x = sin(om_x*x);
  double si_y = sin(om_y*y);

  H->M[3] += D*si_x*si_y;

  dH1->M[3] += om_x*D*cs_x*si_y;
  d2H1->M[3] += -om_x*om_x*D*si_x*si_y;
  dH2->M[3] += om_y*D*si_x*cs_y;
  d2H2->M[3] += -om_y*om_y*D*si_x*si_y;



}

boost::python::list sin_2D_Ham(double x, double y, boost::python::list params_){ 
// sin_2D hamiltonian in diabatic representation

  MATRIX H(2,2);
  MATRIX dH1(2,2);
  MATRIX d2H1(2,2);
  MATRIX dH2(2,2);
  MATRIX d2H2(2,2);


  int sz = boost::python::len(params_);
  vector<double> params(sz,0.0);
  for(int i=0;i<sz;i++){ params[i] = boost::python::extract<double>(params_[i]);  }

  sin_2D_Ham(x,y,&H,&dH1,&d2H1,&dH2,&d2H2,params);

  boost::python::list res;
  res.append(x);
  res.append(H);
  res.append(dH1);
  res.append(d2H1);
  res.append(dH2);
  res.append(d2H2);

 
  return res;

}

}// namespace libmodels
}// liblibra
