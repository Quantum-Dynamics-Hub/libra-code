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
/**
  \file gwp_ints.cpp
  \brief The file implements the integrals needed for Gaussian Wavepacket propagation
    
*/

#include "gwp.h"

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


int check_dimensions(std::string function_name, MATRIX& R1, MATRIX& P1, MATRIX& R2, MATRIX& P2){

  if(R1.num_of_cols!=1){  cout<<"Error in "<<function_name<<": R1 should have only 1 colomn\n"; exit(0);  }
  if(R2.num_of_cols!=1){  cout<<"Error in "<<function_name<<": R2 should have only 1 colomn\n"; exit(0);  }
  if(P1.num_of_cols!=1){  cout<<"Error in "<<function_name<<": P1 should have only 1 colomn\n"; exit(0);  }
  if(P2.num_of_cols!=1){  cout<<"Error in "<<function_name<<": P2 should have only 1 colomn\n"; exit(0);  }

  if(R1.num_of_rows!=P1.num_of_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors R1 (given "<<R1.num_of_rows<<" ) "
        <<"and P1 (given "<<P1.num_of_rows<<" ) do not match\n"; exit(0);     
  }
  if(R1.num_of_rows!=R2.num_of_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors R1 (given "<<R1.num_of_rows<<" ) "
        <<"and R2 (given "<<R2.num_of_rows<<" ) do not match\n"; exit(0);     
  }
  if(P1.num_of_rows!=P2.num_of_rows){ 
    cout<<"Error in "<<function_name<<": The dimensions of vectors P1 (given "<<P1.num_of_rows<<" ) "
        <<"and P2 (given "<<P2.num_of_rows<<" ) do not match\n"; exit(0);     
  }

  // At this point, we are sure that the dimensions of all input vectors are correct
  // so we can return this dimension

  return R1.num_of_rows;

}

complex<double> gwp_overlap(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar){
/**
  This function computes an overlap of two moving Gaussians  <G_1|G_2>, where:

  G_a(r; R_a, P_a, alp_a, gamma_a) = (2*alp_a/pi)^(Ndof/4) * exp(-alp_a*(r-R_a)^2 + i*(P_a/hbar)*(r-R_a) + i*gamma_a/hbar)

  Look derivations at: https://github.com/alexvakimov/Derivatory/blob/master/Gaussian_wavepackets.pdf

  \param[in] R1, R2 Multidimensional vectors containing the components of position of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] P1, P2 Multidimensional vectors containing the components of momenta of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp The Gaussian width factor. Assumed to be the same for all components of both multidimensional Gaussians
  \param[in] hbar The Planck constant in selected units

  The function returns the value of the overlap - a complex number
*/

  int Ndof = check_dimensions("libgwp::gwp_overlap", R1, P1, R2, P2);

  double re = -0.5*alp* ((R2-R1).T() * (R2-R1)).M[0] - (0.125/(alp*hbar*hbar)) * ((P2-P1).T() * (P2-P1)).M[0];
  double im = -(0.5/hbar)* ((P1+P2).T()*(R2-R1)).M[0] + (gamma2 - gamma1)/hbar; 
//  double nrm = pow((2.0*alp/M_PI), 0.25*Ndof); // normalization factor
  
  return exp(complex<double>(re, im));

}


CMATRIX gwp_coupling(MATRIX& R1, MATRIX& P1, double gamma1, 
                     MATRIX& R2, MATRIX& P2, double gamma2, 
                     double alp, double hbar){
/**
  This function computes a derivative coupling between two moving Gaussians  <G_1|d/dr|G_2>, where:

  G_a(r; R_a, P_a, alp_a, gamma_a) = (2*alp_a/pi)^(Ndof/4) * exp(-alp_a*(r-R_a)^2 + i*(P_a/hbar)*(r-R_a) + i*gamma_a/hbar)

  Look derivations at: https://github.com/alexvakimov/Derivatory/blob/master/Gaussian_wavepackets.pdf

  \param[in] R1, R2 Multidimensional vectors containing the components of position of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] P1, P2 Multidimensional vectors containing the components of momenta of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp The Gaussian width factor. Assumed to be the same for all components of both multidimensional Gaussians
  \param[in] hbar The Planck constant in selected units

  The function returns the derivative coupling vector - a complex vector
*/
 
  int Ndof = check_dimensions("libgwp::gwp_coupling", R1, P1, R2, P2);

  // Overlap part 
  MATRIX* re; re = new MATRIX(Ndof,1);  *re = alp*(R2-R1);
  MATRIX* im; im = new MATRIX(Ndof,1);  *im = (0.5/hbar)*(P1+P2);
  CMATRIX* res; res = new CMATRIX(*re, *im);
  
  complex<double> ovlp  = gwp_overlap(R1, P1, gamma1, R2, P2, gamma2, alp, hbar);
  *res *= ovlp;

  delete re; delete im;
  
  return *res;

}


complex<double> gwp_kinetic(MATRIX& R1, MATRIX& P1, double gamma1, 
                            MATRIX& R2, MATRIX& P2, double gamma2, 
                            double alp, double hbar){
/**
  This function computes the kinetic energy matrix element between two moving Gaussians  <G_1|d^2/dr^2|G_2>, where:

  G_a(r; R_a, P_a, alp_a, gamma_a) = (2*alp_a/pi)^(Ndof/4) * exp(-alp_a*(r-R_a)^2 + i*(P_a/hbar)*(r-R_a) + i*gamma_a/hbar)

  Look derivations at: https://github.com/alexvakimov/Derivatory/blob/master/Gaussian_wavepackets.pdf

  \param[in] R1, R2 Multidimensional vectors containing the components of position of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] P1, P2 Multidimensional vectors containing the components of momenta of the multidimensional Gaussians. Assumed
  to be Ndof x 1 vectors
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp The Gaussian width factor. Assumed to be the same for all components of both multidimensional Gaussians
  \param[in] hbar The Planck constant in selected units

  The function returns the matrix element value - a complex number (scalar)
*/
 
  int Ndof = check_dimensions("libgwp::gwp_kinetic", R1, P1, R2, P2);

  // Overlap part 
  double re = -alp + alp*alp*( (R2-R1).T()*(R2-R1) ).M[0] - ( (P2-P1).T()*(P2-P1) ).M[0]/(hbar*hbar)  ;
  double im = alp* ((P1+P2).T() * (R2-R1)).M[0];

  complex<double> res(re, im);  
  complex<double> ovlp  = gwp_overlap(R1, P1, gamma1, R2, P2, gamma2, alp, hbar);

  res *= ovlp;
  
  return res;

}





}// namespace libgwp
}// namespace libdyn
