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
/**
  \file Hamiltonian.cpp
  \brief The file implements the generic Hamiltonian class
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <stdlib.h>
#endif 
#include "Hamiltonian.h"


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libatomistic{


/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


Hamiltonian::Hamiltonian(){ 
/** Default constructor of base (generic Hamiltonian) class
*/
//cout<<"Base Ham. constructor\n";
}

Hamiltonian::~Hamiltonian(){ 
//cout<<"Base Ham. destructor\n"; 
}

void Hamiltonian::set_rep(int rep_){
/**
  Set a wavefunction representation, which affects the Hamiltonian calculations

  \param[in] _rep The representation. Possible options: 0 (for diabatic) and 1 (for adiabatic)
*/

  rep = rep_;
}


void Hamiltonian::set_params(boost::python::list params_){
/**
  Set the parameters of the Hamiltonian
  So far, the model Hamiltonians are implied
  \param[in] params_ The Python list of double-valued parameters to feed into the Hamiltonian
*/

  int sz = boost::python::len(params_);
  vector<double> tmp_params(sz, 0.0);

  // Now copy input params:
  for(int i=0;i<sz; i++){
    tmp_params[i] = boost::python::extract<double>(params_[i]);
  }

  set_params(tmp_params);

}

void Hamiltonian::set_q(vector<double>& q_){
/**
  Update Hamiltonian coordinates (all are the real-valued scalars)

  \param[in] q The vector of real-valued coordinates to be used for Hamiltonian calculations.
  Note: this also sets the status_dia and status_adi variables to zero, impliying the Hamiltonian is not 
  up to date - which is what we want: since the coordinates are changed, the Hamiltonian must be recomputed
  From the prafmatic point of view, if you call this function - expect slower performance.
*/

  q = q_;
  status_dia = 0;
  status_adi = 0;
}

void Hamiltonian::set_q(boost::python::list q_){
/**
  Update Hamiltonian coordinates (all are real-valued scalars - the components of the Python list) - Python-friendly

  \param[in] q The Python list of real-valued coordinates to be used for Hamiltonian calculations.
  Note: this also sets the status_dia and status_adi variables to zero, impliying the Hamiltonian is not 
  up to date - which is what we want: since the coordinates are changed, the Hamiltonian must be recomputed
  From the prafmatic point of view, if you call this function - expect slower performance.
*/

 
  int sz = boost::python::len(q_);
  vector<double> tmp_q(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_q[i] = boost::python::extract<double>(q_[i]);
  }

  set_q(tmp_q);
}

void Hamiltonian::set_v(vector<double>& v_){
/**
  Update Hamiltonian velocities (all are real-valued scalars)

  \param[in] v The vector of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
  Only status_adi is set to 0, so only adiabatic Hamiltonian is recomputed.
  For future: in fact, we only need to update the vibronic Hamiltonian, so we still may save a lot, when adiabatic
  calculations imply electronic structure calculations
*/

  v = v_;
  status_adi = 0;  // only affects adiabatic computations
}

void Hamiltonian::set_v(boost::python::list v_){
/**
  Update Hamiltonian velocities (all are real-valued scalars -the components of Python list) - Python-friendly 

  \param[in] v The vector of real-valued velocities to be used for Hamiltonian calculations.

  The velocities are only needed for vibronic Hamiltonian (adiabatic representation) calculations. 
  Otherwise, they are not used.
  Only status_adi is set to 0, so only adiabatic Hamiltonian is recomputed.
  For future: in fact, we only need to update the vibronic Hamiltonian, so we still may save a lot, when adiabatic
  calculations imply electronic structure calculations
*/

  int sz = boost::python::len(v_);
  vector<double> tmp_v(sz,0.0);

  for(int i=0;i<sz; i++){
    tmp_v[i] = boost::python::extract<double>(v_[i]);
  }

  set_v(tmp_v);

}


void Hamiltonian::compute(){
/**
  Peform actual Hamiltonian computations (is not up to date)

  The computations of either diabatiatic or adiabatic or both Hamiltonians are invoked, depending
  on the representation set up for this Hamiltonian and on the state of the computations of such 
  Hamiltonians (so, if no change of position/velocity has been made since the last computation of 
  given Hamiltonian, no actually computations will be carryied out, not to do usefull work). Also,
  computations of adiabatic Hamiltonians (if adiabatic representation is set up) may call computation
  of the diabatic Hamiltonians, since they may be required. On the contrary, if the diabatic Hamiltonian
  is selected, the adiabatic Hamiltonian is not updated.

  Note, that just updating momenta and positions will not lead to automatic recomputation of the 
  Hamiltonians and derivatives

*/

  if(rep==0){  compute_diabatic();   }
  else if(rep==1){  compute_adiabatic();  }

}



std::complex<double> Hamiltonian::H(int i,int j){
/**
  Return electronic Hamiltonian matrix element

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] i index of electronic state
  \param[in] j index of electronic state

*/


  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_dia->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( ham_adi->get(i,j), 0.0 );
  }

  return res;
}

std::complex<double> Hamiltonian::dHdq(int i,int j,int n){
/**
  Return the derivative of electronic Hamiltonian matrix element w.r.t. nuclear DOF

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] i index of electronic state
  \param[in] j index of electronic state
  \param[in] n index of nuclear DOF w.r.t. which the differentiation is performed

*/


  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_dia[n]->get(i,j), 0.0 );
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>( d1ham_adi[n]->get(i,j), 0.0 );
  }

  return res;
}


std::complex<double> Hamiltonian::D(int i,int j,int n){
/**
  Return the derivative coupling matrix element w.r.t. nuclear DOF

  The returned coupling depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  D = <i|d/dR_n|j> 

  \param[in] i index of electronic state
  \param[in] j index of electronic state
  \param[in] n index of nuclear DOF w.r.t. which the coupling is computed

*/

  std::complex<double> res(0.0,0.0);

  if(rep==0){    // Diabatic Hamiltonian - real, symmetric => Hermitian
    res = std::complex<double>(0.0,0.0);
  }
  else if(rep==1){    // Adiabatic Hamiltonian - real, symmetric => Hermitian
    if(i!=j){  

//      cout<<"in Hamiltonian::D  ... ham_adi = \n"<<*ham_adi<<endl;

      double dE = (ham_adi->get(j,j) - ham_adi->get(i,i) );
      if(fabs(dE)<1e-10){ dE = 1e-10 * (dE>0.0 ? 1.0 : -1.0); }

      res = std::complex<double>( d1ham_adi[n]->get(i,j)/dE, 0.0 );

    }
  }

  return res;
}

std::complex<double> Hamiltonian::nac(int i,int j){
/**
  Return the nonadiabatic coupling matrix element

  The returned coupling depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  nac = sum_n { dR_n/dt * <i|d/dR_n|j> }

  \param[in] traj Index of the trajectory
  \param[in] i index of electronic state
  \param[in] j index of electronic state

*/


  std::complex<double> res(0.0,0.0);

  for(int n=0;n<nnucl;n++){
    res += D(i,j,n) * v[n]; 
  }
  return res;
}

std::complex<double> Hamiltonian::Hvib(int i,int j){
/**
  Return the vibronic Hamiltonian matrix element

  The returned Hamiltonian depends on the selected representation - can be either diabatic or adiabatic.
  This function does not invoke actual computation - it only returns whatever exists in the internal variables.

  \param[in] i index of electronic state
  \param[in] j index of electronic state
*/
 
  const double hbar = 1.0;  // in atomic units

  std::complex<double> ham_ = H(i,j);
  std::complex<double> nac_ = nac(i,j);

  std::complex<double> res(ham_.real(), ham_.imag() - hbar* nac_.real() );

  return res;
}


}// namespace libhamiltonian_generic
}// namespace libatomistic
}// liblibra

