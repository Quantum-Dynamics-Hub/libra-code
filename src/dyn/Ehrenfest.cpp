/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Ehrenfest.cpp
  \brief The file implements all about Ehrenfest calculations
    
*/

#include "Ehrenfest.h"
#include "Dynamics.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace 
namespace libdyn{


void Ehrenfest(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, nHamiltonian& ham, bp::object py_funct, bp::object params, int rep){
 
  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   

  //============== Nuclear propagation ===================
    
       if(rep==0){  p = p + ham.Ehrenfest_forces_dia(C).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(C).real() * 0.5*dt;  }


  q = q + invM*p*dt;
  ham.compute_diabatic(py_funct, bp::object(q), params);
  ham.compute_adiabatic(1);


       if(rep==0){  p = p + ham.Ehrenfest_forces_dia(C).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(C).real() * 0.5*dt;  }

  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   


}

double Ehrenfest_dia(CMATRIX& C, CMATRIX& H, vector<CMATRIX>& dHdR, vector<double>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the diabatic states

  E = <Psi|H_el|Psi> = C.H() * H * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] H Electronic Hamiltonian in the diabatic basis H = <i|H_el|j>
  \param[in] dHdR A vector of derivatives of the Hamiltonian in the diabatic basis w.r.t. all the DOFs:
             dHdR[n].get(i,j) = <i|dH/dR_n|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = (C.H() * H * C).get(0,0).real();

  if(opt>=1){

    int sz = dHdR.size();  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_dia: The size of the output array for forces is not consistent \
      with the size of the input array of the Hamiltonian derivatives\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){
      f[i] = -(C.H() * dHdR[i] * C).get(0,0).real();
    }

  }// opt >=1 

  return energy;
}

double Ehrenfest_dia(CMATRIX* C, CMATRIX* H, vector<CMATRIX*>& dHdR, vector<double*>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the diabatic states

  E = <Psi|H_el|Psi> = C.H() * H * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] H Electronic Hamiltonian in the diabatic basis H = <i|H_el|j>
  \param[in] dHdR A vector of derivatives of the Hamiltonian in the diabatic basis w.r.t. all the DOFs:
             dHdR[n].get(i,j) = <i|dH/dR_n|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = ((*C).H() * (*H) * (*C)).get(0,0).real();

  if(opt>=1){

    int sz = dHdR.size();  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_dia: The size of the output array for forces is not consistent \
      with the size of the input array of the Hamiltonian derivatives\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){
      *f[i] = -((*C).H() * (*dHdR[i]) * (*C)).get(0,0).real();
    }

  }// opt >=1 

  return energy;
}


double Ehrenfest_dia(CMATRIX& C, Hamiltonian& ham, vector<double>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the diabatic states

  E = <Psi|H_el|Psi> = C.H() * H * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] H Electronic Hamiltonian in the diabatic basis H = <i|H_el|j>
  \param[in] dHdR A vector of derivatives of the Hamiltonian in the diabatic basis w.r.t. all the DOFs:
             dHdR[n].get(i,j) = <i|dH/dR_n|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = 0.0; // (C.H() * CMATRIX(ham.get_ham_dia()) * C).get(0,0).real();

  if(opt>=1){

    int sz = ham.nnucl;  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_dia: The size of the output array for forces is not consistent \
      with the number of nuclear DOFs in the input Hamiltonian\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){
      f[i] = 0.0; //-(C.H() * CMATRIX(ham.get_d1ham_dia(i)) * C).get(0,0).real();
    }

  }// opt >=1 

  return energy;
}



double Ehrenfest_adi(CMATRIX& C, CMATRIX& E, vector<CMATRIX>& dEdR, vector<CMATRIX>& D, vector<double>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the adiabatic states

  E = <Psi|H_el|Psi> = C.H() * E * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] E Electronic Hamiltonian in the adiabatic basis
  \param[in] dEdR Derivatives of the adiabatic energies w.r.t. all the DOFs:
             dEdR[n].get(i,i) =  d/dR_n <i|H_el|i>
  \param[in] D A vector of the derivative couplings the Hamiltonian w.r.t. all the DOFs <i|d/dR_alp|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = (C.H() * E * C).get(0,0).real();

  if(opt>=1){

    int sz = D.size();  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_adi: The size of the output array for forces is not consistent \
      with the size of the input array of the derivative couplings\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){
      f[i] = -(C.H() * (dEdR[i] + (D[i]*E-E*D[i])) * C).get(0,0).real();
    }

  }// opt>=1

  return energy;
}


double Ehrenfest_adi(CMATRIX* C, CMATRIX* E, vector<CMATRIX*>& dEdR, vector<CMATRIX*>& D, vector<double*>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the adiabatic states

  E = <Psi|H_el|Psi> = C.H() * E * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] E Electronic Hamiltonian in the adiabatic basis
  \param[in] dEdR Derivatives of the adiabatic energies w.r.t. all the DOFs:
             dEdR[n].get(i,i) =  d/dR_n <i|H_el|i>
  \param[in] D A vector of the derivative couplings the Hamiltonian w.r.t. all the DOFs <i|d/dR_alp|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = ((*C).H() * (*E) * (*C)).get(0,0).real();

  if(opt>=1){

    int sz = D.size();  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_adi: The size of the output array for forces is not consistent \
      with the size of the input array of the derivative couplings\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){
      *f[i] = -( (*C).H() * (*dEdR[i] + ((*D[i]) * (*E) - (*E) * (*D[i]))) * (*C)).get(0,0).real();
    }

  }// opt>=1

  return energy;
}

double Ehrenfest_adi(CMATRIX& C, Hamiltonian& ham, vector<double>& f, int opt){
/**
  \brief Compute Ehrenfest energy and forces of the coherent superposition:

  \Psi(t) = sum_i { C_i(t) * |i> }, where |i> - are the adiabatic states

  E = <Psi|H_el|Psi> = C.H() * E * C

  \param[in] C Time-dependent amplitudes of the basis functions
  \param[in] E Electronic Hamiltonian in the adiabatic basis
  \param[in] dEdR Derivatives of the adiabatic energies w.r.t. all the DOFs:
             dEdR[n].get(i,i) =  d/dR_n <i|H_el|i>
  \param[in] D A vector of the derivative couplings the Hamiltonian w.r.t. all the DOFs <i|d/dR_alp|j>
  \param[out] f A vector of Ehrenfest forces:  f[i] = -dE/dR[i], where E is defined above
  \param[in] opt An option that controls whether to compute Ehrenfest forces (1) or not (0)
*/

  double energy = 0.0; //(C.H() * CMATRIX(ham.get_ham_adi()) * C).get(0,0).real();
/*
  if(opt>=1){

    int sz = D.size();  /// The number of nuclear DOFs
    if(f.size()!=sz){
      cout<<"Error in Ehrenfest_adi: The size of the output array for forces is not consistent \
      with the size of the input array of the derivative couplings\nExiting...";
      exit(0);
    }

    // Do the forces
    for(int i=0;i<sz;i++){

      f[i] = -(C.H() * (dEdR[i] + (D[i]*E-E*D[i])) * C).get(0,0).real();
    }

  }// opt>=1
*/
  return energy;
}



}// namespace libdyn
}// liblibra
