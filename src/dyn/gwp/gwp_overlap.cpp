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
  \file gwp_overlap.cpp
  \brief The file implements the overlap integrals needed for Gaussian Wavepacket propagation
    
*/

#include "gwp.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


complex<double> gwp_overlap(double q1, double p1, double gamma1, double alp1,
                            double q2, double p2, double gamma2, double alp2){
/**
  This function computes an overlap of two moving 1D Gaussians  <G_1|G_2>, where:

  G_a(r; r_a, p_a, alp_a, gamma_a) = (2*alp_a/pi)^(1/4) * exp(-alp_a*(r-r_a)^2 + i*(p_a/hbar)*(r-r_a) + i*gamma_a/hbar)

  Atomic units are assumed, so hbar = 1

  Look derivations at: https://github.com/compchem-cybertraining/derivatory/tree/master/1_gaussian_wavepackets/1_matrix_elements

  \param[in] q1, q2 Coordinates of the Gaussians in a given dimension
  \param[in] p1, p2 Momenta of the Gaussians in a given dimension
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp1, alp2 The Gaussian width factors for given dimension. 


  The function returns the value of the overlap - a complex number
*/

  if(alp1 <= 0.0){  cout<<"Error: alp1 should be positive.\nExiting...\n"; exit(0); }
  if(alp2 <= 0.0){  cout<<"Error: alp2 should be positive.\nExiting...\n"; exit(0); }

  double alp = alp1 + alp2;
  double p_weighted = (alp2 * p1 + alp1 * p2 )/alp; 
  double dq = q2 - q1;
  double dp = p2 - p1;
  double phi = gamma2 - gamma1 - p_weighted * dq;

  double pref1 = alp1 * alp2/alp;

  double res = exp(-pref1 *dq*dq - 0.25*dp*dp/alp );

 
  return pow(4.0*pref1/alp, 0.25) * res * complex<double>(cos(phi), sin(phi));

}


complex<double> gwp_overlap(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2){
/**
  This function computes an overlap of two moving N-dimensional Gaussians  <G_1|G_2>, where:

  G_a(r; r_a, p_a, alp_a, gamma_a) = \product_s^Ndof { (2*alp_a/pi)^(1/4) * exp(-alp_a*(r-r_a)^2 + i*(p_a/hbar)*(r-r_a) + i*gamma_a/hbar)}_s

  Atomic units are assumed, so hbar = 1

  Look derivations at: https://github.com/compchem-cybertraining/derivatory/tree/master/1_gaussian_wavepackets/1_matrix_elements

  \param[in] q1, q2 Coordinates of the Gaussians in a given dimension
  \param[in] p1, p2 Momenta of the Gaussians in a given dimension
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp1, alp2 The Gaussian width factors for given dimension. 

  The function returns the value of the overlap - a complex number
*/

  complex<double> res(1.0, 0.0);


  int Ndof = check_dimensions("libgwp::gwp_overlap", q1, p1, q2, p2);

  for(int  i=0; i<Ndof; i++){
      res *= gwp_overlap(q1.get(i, 0), p1.get(i, 0), gamma1.get(i,0), alp1.get(i,0),
                         q2.get(i, 0), p2.get(i, 0), gamma2.get(i,0), alp2.get(i,0) );
  }

  return res;

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



}// namespace libgwp
}// namespace libdyn
}// liblibra


