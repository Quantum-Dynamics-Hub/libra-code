/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file gwp.cpp
  \brief The file implements various auxiliary functions for Gaussian Wavepacket propagation
    
*/

#include "gwp.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


int check_dimensions(std::string function_name, MATRIX& R1, MATRIX& P1, MATRIX& R2, MATRIX& P2){

  if(R1.n_cols!=1){  cout<<"Error in "<<function_name<<": R1 should have only 1 colomn\n"; exit(0);  }
  if(R2.n_cols!=1){  cout<<"Error in "<<function_name<<": R2 should have only 1 colomn\n"; exit(0);  }
  if(P1.n_cols!=1){  cout<<"Error in "<<function_name<<": P1 should have only 1 colomn\n"; exit(0);  }
  if(P2.n_cols!=1){  cout<<"Error in "<<function_name<<": P2 should have only 1 colomn\n"; exit(0);  }

  if(R1.n_rows!=P1.n_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors R1 (given "<<R1.n_rows<<" ) "
        <<"and P1 (given "<<P1.n_rows<<" ) do not match\n"; exit(0);     
  }
  if(R1.n_rows!=R2.n_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors R1 (given "<<R1.n_rows<<" ) "
        <<"and R2 (given "<<R2.n_rows<<" ) do not match\n"; exit(0);     
  }
  if(P1.n_rows!=P2.n_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors P1 (given "<<P1.n_rows<<" ) "
        <<"and P2 (given "<<P2.n_rows<<" ) do not match\n"; exit(0);     
  }

  // At this point, we are sure that the dimensions of all input vectors are correct
  // so we can return this dimension

  return R1.n_rows;

}


complex<double> gwp_value(MATRIX& r, MATRIX& R, MATRIX& P, double gamma,  double alp, double hbar){
/**
  This function computes the value of the Gaussian at a given point r

  G_a(r; R_a, P_a, alp_a, gamma_a) = (2*alp_a/pi)^(Ndof/4) * exp(-alp_a*(r-R_a)^2 + i*(P_a/hbar)*(r-R_a) + i*gamma_a/hbar)

  Look derivations at: https://github.com/alexvakimov/Derivatory/blob/master/Gaussian_wavepackets.pdf

  \param[in] R Multidimensional vector containing the components of position of the multidimensional Gaussians in all dimensions. Assumed
  to be Ndof x 1 vectors
  \param[in] P Multidimensional vector containing the components of momentum of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] gamma The phase factors of the overall Gaussian
  \param[in] alp The Gaussian width factor. Assumed to be the same for all components of both multidimensional Gaussians
  \param[in] hbar The Planck constant in selected units

  The function returns the value of the Gaussian function - a complex number
*/

  int Ndof = check_dimensions("libgwp::gwp_overlap", r, P, R, P);

  double re = -alp* ((r-R).T() * (r-R)).M[0] ;
  double im =  ( (P.T()*(r-R)).M[0] + gamma)/hbar; 
  double nrm = pow((2.0*alp/M_PI), 0.25*Ndof); // normalization factor
  
  return nrm*exp(complex<double>(re, im));

}



}// namespace libgwp
}// namespace libdyn
}// liblibra


