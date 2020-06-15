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
  \file gwp_kinetic.cpp
  \brief The file implements the kinetic energy integrals needed for Gaussian Wavepacket propagation
    
*/

#include "gwp.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{




complex<double> gwp_kinetic(double q1, double p1, double gamma1, double alp1,
                            double q2, double p2, double gamma2, double alp2){
/**
  This function computes a second derivative matrix element of two moving 1D Gaussians  <G_1|d^2/dr^2|G_2>, where:

  G_a(r; r_a, p_a, alp_a, gamma_a) = (2*alp_a/pi)^(1/4) * exp(-alp_a*(r-r_a)^2 + i*(p_a/hbar)*(r-r_a) + i*gamma_a/hbar)

  Atomic units are assumed, so hbar = 1

  Look derivations at: https://github.com/compchem-cybertraining/derivatory/tree/master/1_gaussian_wavepackets/1_matrix_elements

  \param[in] q1, q2 Coordinates of the Gaussians in a given dimension
  \param[in] p1, p2 Momenta of the Gaussians in a given dimension
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp1, alp2 The Gaussian width factors for given dimension. 


  The function returns the value of the second derivative matrix element - a complex number
*/

  if(alp1 <= 0.0){  cout<<"Error: alp1 should be positive.\nExiting...\n"; exit(0); }
  if(alp2 <= 0.0){  cout<<"Error: alp2 should be positive.\nExiting...\n"; exit(0); }

  double alp = alp1 + alp2;
  double pref1 = 2.0 *alp1 * alp2 / alp; 
  double p_weighted = (alp2 * p1 + alp1 * p2 )/alp; 
  double dq = q2 - q1;
 
  complex<double> ovlp = gwp_overlap(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2);

  return complex<double>( -pref1 + pref1 * pref1 * dq * dq - p_weighted * p_weighted,  2.0 * pref1 * p_weighted * dq) * ovlp;

}


complex<double> gwp_kinetic(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                            MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2){
/**
  This function computes a second derivative matrix elements of two moving N-dimensional Gaussians  <G_1|d^2/dr^2|G_2>, where:

  G_a(r; r_a, p_a, alp_a, gamma_a) = \product_s^Ndof { (2*alp_a/pi)^(1/4) * exp(-alp_a*(r-r_a)^2 + i*(p_a/hbar)*(r-r_a) + i*gamma_a/hbar)}_s

  Atomic units are assumed, so hbar = 1

  Look derivations at: https://github.com/compchem-cybertraining/derivatory/tree/master/1_gaussian_wavepackets/1_matrix_elements

  \param[in] q1, q2 Coordinates of the Gaussians in a given dimension
  \param[in] p1, p2 Momenta of the Gaussians in a given dimension
  \param[in] gamma1, gamma2 The phase factors of the overall Gaussians
  \param[in] alp1, alp2 The Gaussian width factors for given dimension. 

  The function returns the sum of the second derivative matrix elements over all DOFs - one complex number
*/


  int Ndof = check_dimensions("libgwp::gwp_coupling", q1, p1, q2, p2);

  complex<double> res(0.0, 0.0);

  for(int  i=0; i<Ndof; i++){

    if(alp1.get(i,0) <= 0.0){  cout<<"Error: alp1 should be positive.\nExiting...\n"; exit(0); }
    if(alp2.get(i,0) <= 0.0){  cout<<"Error: alp2 should be positive.\nExiting...\n"; exit(0); }


    double alp = alp1.get(i,0) + alp2.get(i,0);
    double pref1 = 2.0 *alp1.get(i,0) * alp2.get(i,0) / alp; 
    double p_weighted = (alp2.get(i,0) * p1.get(i,0) + alp1.get(i,0) * p2.get(i,0) )/alp; 
    double dq = q2.get(i,0) - q1.get(i,0);
 
    res += complex<double>( -pref1 + pref1 * pref1 * dq * dq - p_weighted * p_weighted,  2.0 * pref1 * p_weighted * dq);

  }  

  res *= gwp_overlap(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2);

  return res;

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
  double re = -alp + alp*alp*( (R2-R1).T()*(R2-R1) ).M[0] - 0.25 * ( (P2+P1).T()*(P2+P1) ).M[0]/(hbar*hbar)  ;
  double im = alp* ((P1+P2).T() * (R2-R1)).M[0];

  complex<double> res(re, im);  
  complex<double> ovlp  = gwp_overlap(R1, P1, gamma1, R2, P2, gamma2, alp, hbar);

  res *= ovlp;
  
  return res;

}





}// namespace libgwp
}// namespace libdyn
}// liblibra


