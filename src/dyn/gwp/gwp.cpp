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


double gwp_product_decomposition(double q1, double p1, double gamma1, double alp1,
                                 double q2, double p2, double gamma2, double alp2,
                                 double& q, double& p, double& gamma, double& alp 
                              ){
/**

 This function computes the parameters `q`, `p`, `gamma`, `alp`, and `prefactor` in a Gaussian wavepacket decompositions:

 G(q; q1, p1, gamma1, alp1) *  G(q; q2, p2, gamma2, alp2) =  prefactor * G(q; q, p, gamma, alp)

 This is a version for 1D case

**/


  alp = alp1 + alp2;
  p = p1 + p2;
  q = (alp1 * q1 + alp2 * q2) / alp;
  double dq = q2 - q1;
  gamma = gamma1 + gamma2 + dq * (alp2 * p1 - alp1 * p2) / alp;
  
  double prefactor = pow(2.0 * alp1 * alp2 / alp, 0.25) * exp( - (alp1 * alp2 / alp) * dq*dq ); 
  
  return prefactor;
}


double gwp_product_decomposition(MATRIX& q1, MATRIX& p1, MATRIX& gamma1, MATRIX& alp1,
                                 MATRIX& q2, MATRIX& p2, MATRIX& gamma2, MATRIX& alp2,
                                 MATRIX& q,  MATRIX& p,  MATRIX& gamma,  MATRIX& alp
                              ){
/**

 This function computes the parameters `q`, `p`, `gamma`, `alp`, and `prefactor` in a Gaussian wavepacket decompositions:

 G(q; q1, p1, gamma1, alp1) *  G(q; q2, p2, gamma2, alp2) =  prefactor * G(q; q, p, gamma, alp)

 This is a version for 1D case

**/


  int Ndof = check_dimensions("libgwp::gwp_product_decomposition", q1, p1, q2, p2);
             check_dimensions("libgwp::gwp_product_decomposition", q,  p,  q,  p);


  double prefactor = 1.0;
  double pref1 = 1.0;
  double pref2 = 0.0;

  for(int i=0; i<Ndof; i++){

    alp.set(i, 0, alp1.get(i,0) + alp2.get(i,0));
    p.set(i, 0, p1.get(i,0) + p2.get(i,0));
    q.set(i, 0, (alp1.get(i,0) * q1.get(i,0) + alp2.get(i,0) * q2.get(i,0)) / alp.get(i,0) );
    double dq = q2.get(i,0) - q1.get(i,0);

    gamma.set(i,0, gamma1.get(i,0) + gamma2.get(i,0) + dq * (alp2.get(i,0) * p1.get(i,0) - alp1.get(i,0) * p2.get(i,0)) / alp.get(i,0));      

    double a1a2_over_a = alp1.get(i,0) * alp2.get(i,0) / alp.get(i,0);
    pref1 *= (2.0 * a1a2_over_a );
    pref2 += a1a2_over_a * dq*dq;   
  
  }

  prefactor = pow(pref1, 0.25) * exp( -pref2 );

  return prefactor;
}




}// namespace libgwp
}// namespace libdyn
}// liblibra


