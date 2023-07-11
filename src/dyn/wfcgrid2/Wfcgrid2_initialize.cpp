/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Wfcgrid2_initialize.cpp
  \brief The file implements functions for initializing the wavefunctions on the grids
    
*/

#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libwfcgrid;

/// libwfcgrid namespace
namespace libwfcgrid2{



CMATRIX Gaussian(vector<double>& x, vector<double>& x0, vector<double>& px0, vector<double>& dx, int init_state, int nstates, complex<double> scl){
/**
  \brief Computes a multi-dimensional Gaussian wavepacket at a given point

  \param[out] wfc [CMATRIX(nstates, 1)] Gaussian amplitude at a given state (zero elsewhere)
  \param[in] x is the point of interest (in many dimensions)
  \param[in] x0 position of the center of the Gaussian wavepacket (in many dimensions)
  \param[in] px0 Momentum of the Gaussian wavepacket (in many dimensions)
  \param[in] dx Spread (distribution width) of the spatial component of the Gaussian wavepacket  (in many dimensions)
  \param[in] init_state Index of the occupied electronic state on which the wavepacket is initialized
  \param[in] nstates The number of electronic states for which to initialize the wavefunction

  For each dimension:

  G(x) = [ (1/(2.0*pi*dx^2))^(1/4) ] * exp(-((x-x0)/(2*dx))^2 + i*(x-x0)*px0)

  P(x) = |G(x)|^2 = [ (1/(2.0*pi*dx^2)^2) ] * exp( -(x-x0)^2/(2*dx^2) )

  That is according to: https://en.wikipedia.org/wiki/Normal_distribution,  
  dx - corresponds to standard deviation (in the classical distribution)
  dx^2 - variance (in the classical distribution)

  That is, this wavepacket would correspond to the probability density (classical coordinates)
  that are generated like this:

  x_i = x_ * dx * rnd.normal()


  Other connections:     
  2*a = 1/2*dx^2 =>  a = 1/(2*dx)^2 , where a is such that:  alpha/2 = a + i*b, where a and b are defined in:

  (1) Heller, E. J. Guided Gaussian Wave Packets. Acc. Chem. Res. 2006, 39, 127–134.  
  (2) Akimov, A. V.; Prezhdo, O. V. Formulation of Quantized Hamiltonian Dynamics in Terms of Natural Variables. J. Chem. Phys. 2012, 137, 224115.

*/


  int idof, st;
  int ndof = x.size();

  // Constants
  CMATRIX wfc(nstates,1);
  const complex<double> one(0.0, 1.0);
  double nrm = 1.0;

  for(idof=0; idof<ndof; idof++){
      nrm *= pow((1.0/(2.0*M_PI*dx[idof]*dx[idof])),0.25);
  }

  for(st=0; st<nstates; st++){

      if(st==init_state){

          double c1, c2;
          c1 = 0.0; c2 = 0.0;

          for(idof=0; idof<ndof; idof++){

              double deltx = 0.5*(x[idof] - x0[idof])/dx[idof];

              c1 += -deltx*deltx;
              c2 += px0[idof]*(x[idof] - x0[idof]);
          }

          wfc.M[st] = exp(c1) * (cos(c2)+one*sin(c2));

      }
      else{ wfc.M[st] = complex<double>(0.0, 0.0); }

  }// for st

  wfc *= nrm * scl;

  return wfc;


}



CMATRIX HO(vector<double>& x, vector<double>& x0, vector<double>& px0, 
           vector<double>& alpha, int init_state, int nstates, vector<int>& nu, complex<double> scl){
/**
  \brief Compute n-D Harmonic Oscillator (HO) wavefunction 

  \param[out] wfc [CMATRIX(nstates,1)] Gaussian amplitude at a given state (zero elsewhere)
  \param[in] x is the point of interest (in many dimensions)
  \param[in] x0 Position of the center of the HO basis function in each dimension
  \param[in] px0 Momentum of the HO basis wavepacket in each dimension
  \param[in] alpha The parameter related to the reference HO Hamiltonian:  alpha = sqrt(k*mu/hbar^2)
             Where:  H = -hbar^2 /(2*mu[k]) d^2/dx^2  + 1/2 * k * (x[k]-x0[k])^2 , for each dimension k
  \param[in] init_state Index of the electronic state to which the wavepacket is added
  \param[in] nstates The total number of states
  \param[in] nu Vibrational quantum number of the HO basis function for each dimension
  \param[in] scl The rescaling factor

  The function is:

  weight * prod_{k} [  HO_nu(x[k]-x0[k]) * exp( i*px[k]*(x[k]-x0[k])) ]

  HO_nu(x) = N_nu * H_nu(sqrt(alpha) * (x)) * exp(-1/2 * alpha * x^2 ) 

  N_nu = 1/sqrt(2^nu * nu!)  * (alpha/pi)^(1/4) - normalization factor

  H_nu(ksi):  H_{nu+1}(ksi) - 2 ksi*H_{nu}(ksi) + 2*nu *H_{nu-1}(ksi) = 0  - Hermite polynomial

*/

  int idof, st;
  int ndof = x.size();

  // Constants
  CMATRIX wfc(nstates,1);
  const complex<double> one(0.0, 1.0);
  double nrm = 1.0;

  for(idof=0; idof<ndof; idof++){
      nrm *= (1.0/sqrt(pow(2.0, nu[idof]) * FACTORIAL(nu[idof])) ) * pow((alpha[idof]/M_PI),0.25);
  }


  double H, dH, ksi,ksi2, c2;

  for(st=0; st<nstates; st++){

      if(st==init_state){

        c2 = 0.0;
        ksi2 = 0.0;
        wfc.M[st] = complex<double>(1.0, 0.0);

        for(idof=0; idof<ndof; idof++){
 
          ksi = sqrt(alpha[idof]) * (x[idof] - x0[idof]);
          ksi2 += ksi*ksi;

          HERMITE(nu[idof], ksi, H, dH);

          c2 += px0[idof]*(x[idof] - x0[idof]);

          wfc.M[st] *= H;

        }// for idof

        wfc.M[st] *= exp(-0.5*ksi2) * (cos(c2)+one*sin(c2));

      }// st == init_state
      else{ wfc.M[st] = complex<double>(0.0, 0.0); }

  }// for st

  wfc *= nrm * scl;

  return wfc;


}





void Wfcgrid2::add_wfc_Gau(vector<double>& x0, vector<double>& px0, vector<double>& dx0, int init_state, complex<double> weight, int rep){
/**      
  \brief Add an n-D wavefunction wavepacket to the grid
  \param[in] x0 Position of the center of the Gaussian wavepacket for each dimension
  \param[in] px0 Momentum of the Gaussian wavepacket for each dimension
  \param[in] dx0 Spread (distribution width) of the spatial component of the Gaussian wavepacket for each dimension
  \param[in] init_state Index of the electronic state on which the wavepacket is initialized 
  \param[in] weight - the weight with which the Gaussian will be added to this point
  \param[in] rep - representation ( 0 - diabatic, 1 - adiabatic)

  G(x) = weight * prod_{k=0^ndof} [ (1/(2.0*pi*dx0[k]^2))^(1/4) ] * exp(-((x[k]-x0[k])/(2*dx0[k]))^2 + i*(x[k]-x0[k])*px0[k])
  here, x - in an ndof-dimensional vector
   
*/

  for(int i=0;i<Npts;i++){      
      vector<double> x(ndof, 0.0);

      for(int idof=0; idof<ndof; idof++){
          int ipt = gmap[i][idof];          /// index of the point in that dof
          x[idof] = rgrid[idof]->M[ipt];
      }

      if(rep==0){
        PSI_dia[i] += Gaussian(x, x0, px0, dx0, init_state, nstates, weight);
      }
      else if(rep==1){
        PSI_adi[i] += Gaussian(x, x0, px0, dx0, init_state, nstates, weight);
      }

 
  }

  cout<<"Added a Gaussian to the grid\n";

}// add_wfc_Gau






void Wfcgrid2::add_wfc_HO(vector<double>& x0, vector<double>& px0, vector<double>& alpha, int init_state, vector<int>& nu, complex<double> weight, int rep){
/**      
  \brief Add an n-D HO wavefunction to the grid

  \param[out] wfc [CMATRIX(1, nstates)] Gaussian amplitude at a given state (zero elsewhere)
  \param[in] x is the point of interest (in many dimensions)
  \param[in] x0 Position of the center of the HO basis function in each dimension
  \param[in] px0 Momentum of the HO basis wavepacket in each dimension
  \param[in] alpha The parameter related to the reference HO Hamiltonian:  alpha = sqrt(k*mu/hbar^2)
             Where:  H = -hbar^2 /(2*mu[k]) d^2/dx^2  + 1/2 * k * (x[k]-x0[k])^2 , for each dimension k
  \param[in] init_state Index of the electronic state to which the wavepacket is added
  \param[in] nu Vibrational quantum number of the HO basis function for each dimension
  \param[in] scl The rescaling factor
  \param[in] rep - representation ( 0 - diabatic, 1 - adiabatic)

  The function is:

  weight * prod_{k} [  HO_nu(x[k]-x0[k]) * exp( i*px[k]*(x[k]-x0[k])) ]

  HO_nu(x) = N_nu * H_nu(sqrt(alpha) * (x)) * exp(-1/2 * alpha * x^2 ) 

  N_nu = 1/sqrt(2^nu * nu!)  * (alpha/pi)^(1/4) - normalization factor

  H_nu(ksi):  H_{nu+1}(ksi) - 2 ksi*H_{nu}(ksi) + 2*nu *H_{nu-1}(ksi) = 0  - Hermite polynomial
   
*/

  for(int i=0;i<Npts;i++){      
      vector<double> x(ndof, 0.0);

      for(int idof=0; idof<ndof; idof++){
          int ipt = gmap[i][idof];          /// index of the point in that dof
          x[idof] = rgrid[idof]->M[ipt];
      }

      if(rep==0){
        PSI_dia[i] += HO(x, x0, px0, alpha, init_state, nstates, nu, weight);
      }
      else if(rep==1){
        PSI_adi[i] += HO(x, x0, px0, alpha, init_state, nstates, nu, weight);
      }

 
  }

  cout<<"Added a Harmonic oscillator eigenfunction to the grid\n";

}// add_wfc_HO



void Wfcgrid2::add_wfc_ARB(bp::object py_funct, bp::object params, int rep){
/**
  \brief Initialize a nd-D wavefunction according the external Python function

  \param[in] py_funct - the name of the Python-defined function. Expectations is that the 
  function returns a complex-values result for different states
  \param[in] params - parameters of that function
  \param[in] rep - representation ( 0 - diabatic, 1 - adiabatic)

*/

  CMATRIX res(nstates, 1);

  for(int i=0;i<Npts;i++){      
      MATRIX x(ndof, 1);

      for(int idof=0; idof<ndof; idof++){
          int ipt = gmap[i][idof];          /// index of the point in that dof
          x.set(idof, 0, rgrid[idof]->M[ipt] );
      }

      res = bp::extract< CMATRIX >( py_funct(x, params) );

      if(rep==0){
        PSI_dia[i] += res;
      }
      else if(rep==1){
        PSI_adi[i] += res;
      }
 
  }

  cout<<"Added a custom eigenfunction to the grid\n";

 
}




}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra

