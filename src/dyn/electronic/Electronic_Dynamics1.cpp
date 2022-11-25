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
/**
  \file Electronic_Dynamics1.cpp
  \brief The file implements the methods for solving TD-SE (electronic dynamics)
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <cmath>
#endif 

#include "Electronic.h"
#include "../../math_meigen/libmeigen.h"
#include "../../math_operators/liboperators.h"

/// liblibra namespace
namespace liblibra{


using namespace libmeigen;
using namespace libnhamiltonian;


/// libdyn namespace 
namespace libdyn{

/// libelectronic namespace 
namespace libelectronic{


using liboperators::rotate;


void Electronic::project_out(int i, int renorm_flag){
/**
  \brief Projects the state i out and (optionally) renormalize the wavefunction

  \param[in] i  The index of the state to be projected out
  \param[in] renorm_flag  The parameter to control wavefunction renormalization: 1 - yes, 0 - no

  Warning: "istate" variable is not changed
*/ 

  int j;

  // Project out (annihilate the population of) state i
  q[i] = 0.0;  p[i] = 0.0;

  // Compute the norm after the state projection
  double nrm = 0.0;
  for(j=0;j<nstates;j++){ 
    nrm += (q[j]*q[j] + p[j]*p[j]); 
  } nrm = sqrt(nrm);

  // Optionally, re-normalize the rest of the wavefunction 
  if(renorm_flag==1){
    for(j=0;j<nstates;j++){  q[j] /= nrm; p[j] /= nrm; }
  }
   
}

void Electronic::project_out(int i){
/**
  \brief Projects the state i out and renormalize the wavefunction

  \param[in] i  The index of the state to be projected out

  Warning: "istate" variable is not changed
*/ 

  this->project_out(i,1);
}


void Electronic::collapse(int i, int phase_flag){
/**
  \brief Collapses the wavefunction onto the state i 

  \param[in] i  The index of the state onto which the wfc is collapsed
  \param[in] phase_flag Controls how the phase is to be treated:
   phase_flag = 0 - do not care about phase - destroy it
   phase_flag = 1 - preserve the phase

  This also affects the "istate" variable - it is changed
*/ 

  // Change the state identity
  istate = i;

  // Take care of the amplitudes:
  for(int j=0;j<nstates;j++){
    if(j==i){
      if(phase_flag==0){  q[i] = 1.0; p[i] = 0.0;   }
      else if(phase_flag==1){
        double popi = q[i]*q[i] + p[i]*p[i];  // original population of state i
        double nrm = sqrt(popi);

        if(popi>0.0){ q[i] /= nrm; p[i] /= nrm; }
        else{ q[i] = 1.0; p[i] = 0.0; }
       
      }// phase_flag==1
    }// j==i
    else{  q[j] = p[j] = 0.0; }
  }// for j

}

void Electronic::collapse(int i){
/**
  \brief Collapses the wavefunction onto the state i and preserves the phase

  \param[in] i  The index of the state onto which the wfc is collapsed

  This also affects the "istate" variable - it is changed
*/ 
  this->collapse(i, 1);
}


/**

void Electronic::propagate_electronic(double dt,Hamiltonian* ham){
  \brief Propagate electronic DOF using sequential rotations in the MMTS variables

  Methodologically: This version is based on the Hamiltonian formulation of TD-SE
  This propagator is good for general Hamiltonian - diabatic of adiabatic
  iL = iL_qp + iL_qq + iL_pp
  iL_qp = sum_i,j {}

  API: The member function of the Electronic function, with the pointer to Hamiltonian object

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in] ham The pointer to Hamiltonian object, that affects the dynamics

  libdyn::libelectronic::propagate_electronic(dt,this,ham);
}
*/ 

/**

void Electronic::propagate_electronic(double dt,Hamiltonian& ham){
  \brief Propagate electronic DOF using sequential rotations in the MMTS variables

  Methodologically: This version is based on the Hamiltonian formulation of TD-SE
  This propagator is good for general Hamiltonian - diabatic of adiabatic
  iL = iL_qp + iL_qq + iL_pp
  iL_qp = sum_i,j {}

  API: The member function of the Electronic function, with the reference to Hamiltonian object

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in] ham The reference to Hamiltonian object, that affects the dynamics

  libdyn::libelectronic::propagate_electronic(dt,this,&ham);
}
*/ 

/**

void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham){
  \brief Propagate electronic DOF using sequential rotations in the MMTS variables

  Methodologically: This version is based on the Hamiltonian formulation of TD-SE
  This propagator is good for general Hamiltonian - diabatic of adiabatic
  iL = iL_qp + iL_qq + iL_pp
  iL_qp = sum_i,j {}

  API: A free function that takes Electronic object as the input and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] el The pointer to the Electronic object containing the electronic DOF
  \param[in] ham The pointer to Hamiltonian object, that affects the dynamics


  int i,j;

  double dt_half = 0.5*dt;

  //------------- Phase evolution (adiabatic) ---------------- 
  // exp(iL_qp * dt/2)
  for(i=0;i<el->nstates;i++){
    for(j=0;j<el->nstates;j++){

      rotate(el->p[j],el->q[i], dt_half*ham->Hvib(i,j).real());

    }// for j
  }// for i

  //------------- Population transfer (adiabatic) ----------------
  // exp((iL_qq + iL_pp) * dt/2)
  for(i=0;i<el->nstates;i++){
    for(j=i+1;j<el->nstates;j++){

      rotate(el->q[j],el->q[i], dt_half*ham->Hvib(i,j).imag());
      rotate(el->p[j],el->p[i], dt_half*ham->Hvib(i,j).imag());

    }// for j
  }// for i

  // exp((iL_qq + iL_pp) * dt/2)
  for(i=el->nstates-1;i>=0;i--){
    for(j=el->nstates-1;j>i;j--){


      rotate(el->q[j],el->q[i], dt_half*ham->Hvib(i,j).imag());
      rotate(el->p[j],el->p[i], dt_half*ham->Hvib(i,j).imag());

    }// for j
  }// for i


  //------------- Phase evolution (adiabatic) ----------------
  // exp(iL_qp * dt/2)
  for(i=el->nstates-1;i>=0;i--){
    for(j=el->nstates-1;j>=0;j--){

      rotate(el->p[j],el->q[i], dt_half*ham->Hvib(i,j).real());

    }// for j
  }// for i


}// propagate_electronic

*/ 


void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib){
/**
  \brief Propagate electronic DOF using sequential rotations in the MMTS variables

  Methodologically: This version is based on the Hamiltonian formulation of TD-SE
  This propagator is good for general Hamiltonian - diabatic of adiabatic
  iL = iL_qp + iL_qq + iL_pp
  iL_qp = sum_i,j {}

  API: A free function that takes Electronic object as the input and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] el The reference to the Electronic object containing the electronic DOF
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) 

  This is the Python-friendly function
*/ 


  int i,j;

  double dt_half = 0.5*dt;

  //------------- Phase evolution (adiabatic) ---------------- 
  // exp(iL_qp * dt/2)
  for(i=0;i<el.nstates;i++){
    for(j=0;j<el.nstates;j++){

      rotate(el.p[j],el.q[i], dt_half*Hvib.get(i,j).real());

    }// for j
  }// for i

  //------------- Population transfer (adiabatic) ----------------
  // exp((iL_qq + iL_pp) * dt/2)
  for(i=0;i<el.nstates;i++){
    for(j=i+1;j<el.nstates;j++){

      rotate(el.q[j],el.q[i], dt_half*Hvib.get(i,j).imag());
      rotate(el.p[j],el.p[i], dt_half*Hvib.get(i,j).imag());

    }// for j
  }// for i

  // exp((iL_qq + iL_pp) * dt/2)
  for(i=el.nstates-1;i>=0;i--){
    for(j=el.nstates-1;j>i;j--){


      rotate(el.q[j],el.q[i], dt_half*Hvib.get(i,j).imag());
      rotate(el.p[j],el.p[i], dt_half*Hvib.get(i,j).imag());

    }// for j
  }// for i


  //------------- Phase evolution (adiabatic) ----------------
  // exp(iL_qp * dt/2)
  for(i=el.nstates-1;i>=0;i--){
    for(j=el.nstates-1;j>=0;j--){

      rotate(el.p[j],el.q[i], dt_half*Hvib.get(i,j).real());

    }// for j
  }// for i


}// propagate_electronic



void propagate_electronic_rot(double dt, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  \brief Propagate electronic DOF using sequential rotations in the MMTS variables

  Methodologically: This version is based on the Hamiltonian formulation of TD-SE
  This propagator is good for general Hamiltonian - diabatic of adiabatic
  iL = iL_qp + iL_qq + iL_pp
  iL_qp = sum_i,j {}

  API: A free function that takes Electronic object as the input and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff(nst, nst) The reference to the Electronic object containing the electronic DOF
  \param[in] Hvib The reference to the vibronic Hamiltonian matrix

  This is the Python-friendly function
*/ 


  int i,j;
  double qi,pi, qj, pj;

  double dt_half = 0.5*dt;

  if(Coeff.n_cols!=1){
    cout<<"ERROR in propagate_electronic_rot: The number of columns in the input amplitudes vector\
          should be 1, but "<<Coeff.n_cols<<" is given \n";
    exit(0);
  }
  if(Hvib.n_cols!=Hvib.n_rows){
    cout<<"ERROR in propagate_electronic_rot: The number of columns and rows of the input Hamiltonian\
          should be equal to each other, but the numbers are: \n";
    cout<<"num_of_rows = "<<Hvib.n_rows<<"\n";
    cout<<"num_of_cols = "<<Hvib.n_cols<<"\n";
    exit(0);
  }
  if(Hvib.n_cols!=Coeff.n_rows){
    cout<<"ERROR in propagate_electronic_rot: The number of columns of the input Hamiltonian\
          should be equal to the number of rows of the amplitudes vector, but what is given: \n";
    cout<<"Coeff.num_of_rows = "<<Coeff.n_rows<<"\n";
    cout<<"Hvib.num_of_cols = "<<Hvib.n_cols<<"\n";
    exit(0);
  }


  int nstates = Coeff.n_rows;

  //------------- Phase evolution (adiabatic) ---------------- 
  // exp(iL_qp * dt/2)
  for(i=0;i<nstates;i++){
    for(j=0;j<nstates;j++){

      qi = Coeff.get(i,0).real();  pi = Coeff.get(i,0).imag();
      qj = Coeff.get(j,0).real();  pj = Coeff.get(j,0).imag();

      rotate(pj, qi, dt_half*Hvib.get(i,j).real());

      Coeff.set(i,0, complex<double>(qi, pi));
      Coeff.set(j,0, complex<double>(qj, pj));

    }// for j
  }// for i

  //------------- Population transfer (adiabatic) ----------------
  // exp((iL_qq + iL_pp) * dt/2)
  for(i=0;i<nstates;i++){
    for(j=i+1;j<nstates;j++){

      qi = Coeff.get(i,0).real();  pi = Coeff.get(i,0).imag();
      qj = Coeff.get(j,0).real();  pj = Coeff.get(j,0).imag();

      rotate(qj, qi, dt_half*Hvib.get(i,j).imag());
      rotate(pj, pi, dt_half*Hvib.get(i,j).imag());

      Coeff.set(i,0, complex<double>(qi, pi));
      Coeff.set(j,0, complex<double>(qj, pj));


    }// for j
  }// for i

  // exp((iL_qq + iL_pp) * dt/2)
  for(i=nstates-1;i>=0;i--){
    for(j=nstates-1;j>i;j--){

      qi = Coeff.get(i,0).real();  pi = Coeff.get(i,0).imag();
      qj = Coeff.get(j,0).real();  pj = Coeff.get(j,0).imag();

      rotate(qj, qi, dt_half*Hvib.get(i,j).imag());
      rotate(pj, pi, dt_half*Hvib.get(i,j).imag());

      Coeff.set(i,0, complex<double>(qi, pi));
      Coeff.set(j,0, complex<double>(qj, pj));

    }// for j
  }// for i


  //------------- Phase evolution (adiabatic) ----------------
  // exp(iL_qp * dt/2)
  for(i=nstates-1;i>=0;i--){
    for(j=nstates-1;j>=0;j--){

      qi = Coeff.get(i,0).real();  pi = Coeff.get(i,0).imag();
      qj = Coeff.get(j,0).real();  pj = Coeff.get(j,0).imag();

      rotate(pj, qi, dt_half*Hvib.get(i,j).real());

      Coeff.set(i,0, complex<double>(qi, pi));
      Coeff.set(j,0, complex<double>(qj, pj));

    }// for j
  }// for i


}// propagate_electronic_rot




void propagate_electronic_eig(double dt, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  Solves the time-dependent Schrodinger equation:

  i*hbar*dc/dt = Hvib*c 
   
  API: A free function that takes electronic DOF in the form of matrix-colomun and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff The reference to the CMATRIX object containing the electronic DOF (coefficient)
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
             it is also assumed to Hermitian - pay attention to how it is constructed!
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued time-dependent matrix, CMATRIX)

  This integrator IS unitary
  This is the Python-friendly function
*/ 


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = Hvib.n_cols;  

  // Transform the Hamiltonian accordingly:
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity();
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0); // eigenvectors
  CMATRIX* Heig; Heig = new CMATRIX(sz, sz);  *Heig = complex<double>(0.0,0.0); // eigenvalues
  CMATRIX* expH;   expH = new CMATRIX(sz, sz);    *expH = complex<double>(0.0,0.0);   


  // Compute the exponential  exp(-i*Hvib*dt)  
  libmeigen::solve_eigen(Hvib, *I, *Heig, *C, 0);  // Hvib_eff * C = I * C * Heig  ==>  Hvib = C * Heig * C.H()


  // Diagonal form expH
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Heig->get(i,i)*dt );
    expH->set(i,i,val);
  }

  // Transform back to the original basis:
  *expH = (*C) * (*expH) * ((*C).H());

  // Propagation
  Coeff = *expH * Coeff;

  
  // Clean temporary memory
  delete expH;  delete Heig;  delete C;  delete I;



}// propagate_electronic

void propagate_electronic(double dt, CMATRIX& Coeff, CMATRIX& Hvib){

  int i,j;
  // Faster and more robust solver:
  int ntraj = Coeff.n_cols;
  int nel = Coeff.n_rows;  
  Electronic el(nel);

  for(j=0; j<ntraj; j++){

    for(i=0; i<nel; i++){ 
      el.q[i] = Coeff.get(i,j).real(); 
      el.p[i] = Coeff.get(i,j).imag();       
    }
    propagate_electronic(dt, el, Hvib);

    for(i=0; i<nel; i++){   
      Coeff.set(i,j, complex<double>(el.q[i], el.p[i]));    
    }

  }

  // Older Default 
  //propagate_electronic_eig(dt, Coeff, Hvib);
}



void propagate_electronic_nonHermitian(double dt, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  Solves the time-dependent Schrodinger equation:

  i*hbar*dc/dt = Hvib*c  with Hvib being non-Hermitian
   
  API: A free function that takes electronic DOF in the form of matrix-colomun and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff The reference to the CMATRIX object containing the electronic DOF (coefficient)
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
             it is also assumed to Hermitian - pay attention to how it is constructed!
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued time-dependent matrix, CMATRIX)

  This is the Python-friendly function

  http://www.pha.jhu.edu/~neufeld/numerical/lecturenotes10.pdf

  A * R = R * D                               
  L * A  = D * L <=> A.H() * L.H() = L.H() * D

  Then  L.H() * R = I 

  So:  D = L.H() * A * R  and   A = R * D * L

*/ 


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = Hvib.n_cols;  

  // Transform the Hamiltonian accordingly:
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity();
  CMATRIX* R; R = new CMATRIX(sz, sz);  *R = complex<double>(0.0, 0.0); // right eigenvectors:  A * R = R * D
  CMATRIX* L; L = new CMATRIX(sz, sz);  *L = complex<double>(0.0, 0.0); // left eigenvectors:   L * A  = D * L <=> A.H() * L.H() = L.H() * D
  CMATRIX* Heig; Heig = new CMATRIX(sz, sz);  *Heig = complex<double>(0.0,0.0); // eigenvalues
  CMATRIX* HvibH; HvibH = new CMATRIX(sz,sz); *HvibH = Hvib.H();
  CMATRIX* expH;   expH = new CMATRIX(sz, sz);    *expH = complex<double>(0.0,0.0);   


  // Solve the right eigenvalue problem
  libmeigen::solve_eigen(Hvib, *I, *Heig, *R, 0);  // Hvib * R = R * Heig  

  // Solve the left eigenvalue problem
  libmeigen::solve_eigen(*HvibH, *I, *Heig, *L, 0);  // Hvib.H() * L.H() = L.H() * Heig  
  *L = (*L).H();


  // Diagonal form expH
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Heig->get(i,i)*dt );
    expH->set(i,i,val);
  }

  // Transform back to the original basis:
  *expH = (*R) * (*expH) * (*L);

  // Propagation
  Coeff = (*expH) * Coeff;

  
  // Clean temporary memory
  delete expH;  delete Heig; delete HvibH; delete R; delete L;  delete I;



}// propagate_electronic






void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, MATRIX& S){
/**
  \brief Propagate electronic DOF in the MMTS variables

  Methodologically: solve i*hbar*S*dc/dt = Hvib*c directly, using Lowdin-type transformation and 
  matrix exponentiation. This solver is good only for time-independent S matrix

  API: A free function that takes Electronic object as the input and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] el The reference to the Electronic object containing the electronic DOF
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
  \param[in] S The reference to the overlap matrix (assumed to be a real-valued matrix, MATRIX)

  This integrator is fully unitary (the norm is conserved exactly)
  This is the Python-friendly function
*/ 

  int i,j;

  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;
  
  MATRIX* I; I = new MATRIX(sz, sz);  I->Init_Unit_Matrix(1.0);
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);

  // Transformation to adiabatic basis
  libmeigen::solve_eigen(S, *I, *Seig, *C, 0);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


  CMATRIX* S_i_half; S_i_half = new CMATRIX(sz, sz);  *S_i_half = complex<double>(0.0,0.0);  // S^{-1/2}
  CMATRIX* S_half;   S_half = new CMATRIX(sz, sz);    *S_half = complex<double>(0.0,0.0);    // S^{1/2}


  // Diagonal form of S^{-1/2} and S^{1/2} matrices
  for(i=0;i<sz;i++){
    complex<double> val = std::sqrt(Seig->get(i,i));

    S_i_half->M[i*sz+i] = 1.0/val;
    S_half->M[i*sz+i] = val;
  }

  // Convert to original basis
  *S_i_half = (*C) * (*S_i_half) * ((*C).H());
  *S_half = (*C) * (*S_half) * ((*C).H());




  // Transform the Hamiltonian accordingly:
  CMATRIX* Hvib_eff;  Hvib_eff = new CMATRIX(sz,sz);
  CMATRIX* coeff;     coeff = new CMATRIX(sz,1);
  for(i=0;i<sz;i++){  coeff->M[i] = complex<double>(el.q[i], el.p[i]);  }


  // Transform Hvib and coefficients
  *Hvib_eff = (*S_i_half) * (Hvib) * (*S_i_half);
  *coeff = (*S_half) * (*coeff); // * (*S_half);      // now these are effective coefficients
  
  // Set up the object
  Electronic el_eff; el_eff = el;
  for(i=0; i<sz; i++){  el_eff.q[i] = coeff->get(i,0).real(); el_eff.p[i] = coeff->get(i,0).imag(); }


  // Now do the standard propagation
  propagate_electronic(dt, el_eff, *Hvib_eff);


  for(i=0;i<sz;i++){  coeff->M[i] = complex<double>(el_eff.q[i], el_eff.p[i]);  }


  // Transform the coefficients back to the original representation:
  *coeff = (*S_i_half) * (*coeff); // * (*S_i_half);      // now these are effective coefficients
  

  // Update "normal" electronic variables
  for(i=0; i<sz; i++){  el.q[i] = coeff->get(i,0).real(); el.p[i] = coeff->get(i,0).imag(); }  


  // Clean temporary memory
  delete coeff;
  delete Hvib_eff;
  delete S_i_half;
  delete S_half;
  delete Seig;
  delete C;
  delete I;



}// propagate_electronic

/**

void Electronic::propagate_electronic(double dt,Hamiltonian& ham, CMATRIX& S){
  \brief Propagate electronic DOF 

  Same as 
  void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S)

  CMATRIX* Hvib; Hvib = new CMATRIX(nstates, nstates);

  for(int i=0;i<nstates;i++){
    for(int j=0;j<nstates;j++){
      Hvib->set(i,j, ham.Hvib(i,j) );
    }
  }
  
  libdyn::libelectronic::propagate_electronic(dt, *this, *Hvib, S);

  delete Hvib;

}

*/ 


void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S){
/**
  \brief Propagate electronic DOF in the MMTS variables

  Methodologically: solve i*hbar*S*dc/dt = Hvib*c directly, using Lowdin-type transformation and 
  matrix exponentiation. This solver is good only for time-independent S matrix

  API: A free function that takes Electronic object as the input and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] el The reference to the Electronic object containing the electronic DOF
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued matrix, CMATRIX)

  This integrator is fully unitary (the norm is conserved exactly)
  This is the Python-friendly function
*/ 

  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Transformation to adiabatic basis
  libmeigen::solve_eigen(S, *I, *Seig, *C, 0);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


  // Diagonal form of S^{-1/2} and S^{1/2} matrices
  CMATRIX* S_i_half; S_i_half = new CMATRIX(sz, sz);  *S_i_half = complex<double>(0.0,0.0);  // S^{-1/2}
  CMATRIX* S_half;   S_half = new CMATRIX(sz, sz);    *S_half = complex<double>(0.0,0.0);    // S^{1/2}
  for(i=0;i<sz;i++){
    complex<double> val = std::sqrt(Seig->get(i,i));

    S_i_half->M[i*sz+i] = 1.0/val;
    S_half->M[i*sz+i] = val;
  }


  // Convert to original basis
  *S_i_half = (*C) * (*S_i_half) * ((*C).H());
  *S_half = (*C) * (*S_half) * ((*C).H());


  // Transform the Hamiltonian accordingly:
  CMATRIX* Hvib_eff;  Hvib_eff = new CMATRIX(sz,sz);
  CMATRIX* coeff;     coeff = new CMATRIX(sz,1);
  for(i=0;i<sz;i++){  coeff->M[i] = complex<double>(el.q[i], el.p[i]);  }


  // Transform Hvib and coefficients
  *Hvib_eff = (*S_i_half) * (Hvib) * (*S_i_half);
  *coeff = (*S_half) * (*coeff);                 // now these are effective coefficients

  
  // Set up the object
  Electronic el_eff; el_eff = el;
  for(i=0; i<sz; i++){  el_eff.q[i] = coeff->get(i,0).real(); el_eff.p[i] = coeff->get(i,0).imag(); }


  // Compute the exponential  exp(-i*Hvib*dt)  
  libmeigen::solve_eigen(*Hvib_eff, *I, *Seig, *C, 0);  // Hvib_eff * C = I * C * Seig  ==>  Hvib = C * Seig * C.H()


  // Diagonal form expH
  CMATRIX* expH;   expH = new CMATRIX(sz, sz);    *expH = complex<double>(0.0,0.0);    // 
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Seig->get(i,i)*dt );
    expH->set(i,i,val);
  }


  // Transform back to the original basis:
  *expH = (*C) * (*expH) * ((*C).H());


  // Propagation
   *coeff = *expH * *coeff;


  // Transform the coefficients back to the original representation:
  *coeff = (*S_i_half) * (*coeff); // * (*S_i_half);      // now these are effective coefficients

  
  // Update "normal" electronic variables
  for(i=0; i<sz; i++){  el.q[i] = coeff->get(i,0).real(); el.p[i] = coeff->get(i,0).imag(); }  


  // Clean temporary memory
  delete expH;
  delete coeff;
  delete Hvib_eff;
  delete S_i_half;
  delete S_half;
  delete Seig;
  delete C;
  delete I;



}// propagate_electronic




void propagate_electronic(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S){
/**
  Solves the generalized time-dependent Schrodinger equation:

  i*hbar*S*dc/dt = Hvib*c 

  Although the vibronic Hamiltonian, Hvib_tilda, computed in the non-orthonormal basis is not Hermitian,
  the theory eventually leads to a Hermitian version of the vibronic Hamiltonian, so we take the Hermitian one
  as the input. See my paper on the details.
 
  API: A free function that takes electronic DOF in the form of matrix-colomun and modifies it

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff The reference to the CMATRIX object containing the electronic DOF (coefficient)
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
             it is also assumed to Hermitian - pay attention to how it is constructed!
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued time-dependent matrix, CMATRIX)

  This integrator IS unitary
  This is the Python-friendly function
*/ 


  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Transformation to adiabatic basis
  libmeigen::solve_eigen(S, *I, *Seig, *C, 0);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


  // Diagonal form of S^{-1/2} and S^{1/2} matrices
  CMATRIX* S_i_half; S_i_half = new CMATRIX(sz, sz);  *S_i_half = complex<double>(0.0,0.0);  // S^{-1/2}
  CMATRIX* S_half;   S_half = new CMATRIX(sz, sz);    *S_half = complex<double>(0.0,0.0);    // S^{1/2}
  for(i=0;i<sz;i++){
    complex<double> val = std::sqrt(Seig->get(i,i));

    S_i_half->M[i*sz+i] = 1.0/val;
    S_half->M[i*sz+i] = val;
  }


  // Convert to original basis
  *S_i_half = (*C) * (*S_i_half) * ((*C).H());
  *S_half = (*C) * (*S_half) * ((*C).H());


  // Transform the Hamiltonian accordingly:
  CMATRIX* Hvib_eff;  Hvib_eff = new CMATRIX(sz,sz);
  CMATRIX* coeff;     coeff = new CMATRIX(sz,1);

  *Hvib_eff = (*S_i_half) * Hvib * (*S_i_half);  // Hermitian part
  *coeff = (*S_half) * Coeff;                    // now these are the effective coefficients

  
  // Compute the exponential  exp(-i*Hvib*dt)  
  libmeigen::solve_eigen(*Hvib_eff, *I, *Seig, *C, 0);  // Hvib_eff * C = I * C * Seig  ==>  Hvib = C * Seig * C.H()


  // Diagonal form expH
  CMATRIX* expH;   expH = new CMATRIX(sz, sz);    *expH = complex<double>(0.0,0.0);    // 
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Seig->get(i,i)*dt );
    expH->set(i,i,val);
  }

  // Transform back to the original basis:
  *expH = (*C) * (*expH) * ((*C).H());

  // Propagation
   *coeff = *expH * *coeff;

  // Transform the coefficients back to the original representation:
  *coeff = (*S_i_half) * (*coeff); // * (*S_i_half);      // convert the effective coefficients back to original representation

  for(i=0;i<Coeff.n_elts;i++){  Coeff.M[i] = coeff->M[i]; } // so we don't allocate new memory for Coeff!!!

  
  // Clean temporary memory
  delete expH;
  delete coeff;
  delete Hvib_eff;
  delete S_i_half;
  delete S_half;
  delete Seig;
  delete C;
  delete I;



}// propagate_electronic




void propagate_electronic_qtag(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S){
/**
  Solves the generalized time-dependent Schrodinger equation:

  i*hbar*S*dc/dt = Hvib*c 
 
  using the approach described in the QTAG paper:


  U = Z exp(-i/hbar *E *dt) * Z^+  * S

  This is a potentially more robust integrator that the other one, where we compute S^+/-{1/2}

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff The reference to the CMATRIX object containing the electronic DOF (coefficient)
  \param[in] ham The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
             it is also assumed to Hermitian - pay attention to how it is constructed!
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued time-dependent matrix, CMATRIX)

*/ 


  int i;
 
  // Let us first diagonalize the overlap matrix Hvib
  int sz = Hvib.n_cols;  
  CMATRIX Z(sz, sz);
  CMATRIX E(sz, sz);

  // Transformation to adiabatic basis
  libmeigen::solve_eigen(Hvib, S, E, Z, 0);  // Hvib * Z = S * Z * E


  // Diagonal form expH
  CMATRIX expE(sz, sz); 
  for(i=0;i<sz;i++){
    double argg = -E.get(i,i).real()*dt;
    double _cs = cos(argg);
    double _si = sin(argg);
    expE.set(i,i, complex<double>(_cs, _si)  );
  }


  // This is a supposedly more efficient form of Coeff = (Z * expE * Z.H() * S) * Coeff
  Coeff = S * Coeff;
  Coeff = Z.H() * Coeff;
  Coeff = expE * Coeff;
  Coeff = Z * Coeff;


}// propagate_electronic_qtag



void propagate_electronic(double dt, CMATRIX& C, nHamiltonian& ham, int rep){

 if(rep==0){  // diabatic

    CMATRIX Hvib(ham.ndia, ham.ndia);  Hvib = ham.get_hvib_dia();
    CMATRIX Sdia(ham.ndia, ham.ndia);  Sdia = ham.get_ovlp_dia();
//    Hvib = Hvib 

    propagate_electronic(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    //propagate_electronic_qtag(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs

    //propagate_electronic_nonHermitian(dt, C, Hvib);

  }

  else if(rep==1){  // adiabatic

    CMATRIX Hvib(ham.nadi, ham.nadi);  Hvib = ham.get_hvib_adi(); 

    propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs
  }

}



void propagate_electronic(double dt, CMATRIX& C, nHamiltonian* ham, int rep){

  if(rep==0){  // diabatic

    CMATRIX Hvib(ham->ndia, ham->ndia);  Hvib = ham->get_hvib_dia();
    CMATRIX Sdia(ham->ndia, ham->ndia);  Sdia = ham->get_ovlp_dia();

    propagate_electronic(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    //propagate_electronic_qtag(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs

  }

  else if(rep==1){  // adiabatic

    CMATRIX Hvib(ham->nadi, ham->nadi);  Hvib = ham->get_hvib_adi(); 
    propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs

  }

}


void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep){

  if(C.n_cols!=ham.size()){
    cout<<"ERROR in void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep): \n";
    cout<<"C.n_cols = "<<C.n_cols<<" is not equal to ham.size() = "<<ham.size()<<"\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  int nst = C.n_rows;
  int ntraj = C.n_cols;
  
  CMATRIX ctmp(nst, 1);

  for(int traj=0; traj<ntraj; traj++){
    ctmp = C.col(traj);
    propagate_electronic(dt, ctmp, ham[traj], rep);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }

  }

}

void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep, int isNBRA){
/**
  if isNBRA == 1, then only the first element of ham (ham[0]) will be used. However, we assume that the C matrix still
  has the full dimensionality, that is it consist of ntraj columns
*/

  if(!isNBRA){  
    if(C.n_cols!=ham.size()){
      cout<<"ERROR in void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep): \n";
      cout<<"C.n_cols = "<<C.n_cols<<" is not equal to ham.size() = "<<ham.size()<<"\n";
      cout<<"Exiting...\n";
      exit(0);
    }
  }

  int nst = C.n_rows;
  int ntraj = C.n_cols;
  
  CMATRIX ctmp(nst, 1);

  for(int traj=0; traj<ntraj; traj++){

    ctmp = C.col(traj);

    int traj1 = traj; 
    if(isNBRA==1){ traj1 = 0; }

    propagate_electronic(dt, ctmp, ham[traj1], rep);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }

  }

}






void propagate_electronic(double dt, CMATRIX& C, nHamiltonian& ham, int rep, int level){

  vector<nHamiltonian*> branches; 
  branches = ham.get_branches(level);

  if(C.n_cols!=branches.size()){
    cout<<"ERROR in void propagate_electronic(double dt, CMATRIX& C, nHamiltonian& ham, int rep, int level): \n";
    cout<<"C.n_cols = "<<C.n_cols<<" is not equal to ham.size() = "<<branches.size()<<"\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  int nst = C.n_rows;
  int ntraj = C.n_cols;
  
  CMATRIX ctmp(nst, 1);

  for(int traj=0; traj<ntraj; traj++){
    ctmp = C.col(traj);
    propagate_electronic(dt, ctmp, branches[traj], rep);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }

  }

}




void grid_propagator(double dt, CMATRIX& Hvib, CMATRIX& S, CMATRIX& U){
/**
  A super-fast version specifically taylored to grid integration (since grid is known,
  the matrices S, Hvib and all the dependent matrices are fixed, so can be computed only once and then
  passed here)

  Solves the generalized time-dependent Schrodinger equation:

  i*hbar*S*dc/dt = Hvib*c 
 
  \param[in] dt The integration time step (also the duration of propagation)
  \param[in] Hvib The reference to the vibronic Hamiltonian matrix (not the Hamiltonian object!) - the complex-valued matrix, CMATRIX
             it is also assumed to be Hermitian - pay attention to how it is constructed!
  \param[in] S The reference to the overlap matrix (assumed to be a complex-valued time-dependent matrix, CMATRIX)

  This is the Python-friendly function
  This function returns a propagator
*/ 

  int sz = S.n_cols;  
  CMATRIX S_half(sz,sz);
  CMATRIX S_i_half(sz,sz);

  // Diagonal form of S^{-1/2} and S^{1/2} matrices
  sqrt_matrix(S, S_half, S_i_half);

  // Transform the Hamiltonian accordingly:
  CMATRIX Hvib_eff(sz,sz);  
  Hvib_eff = S_i_half * Hvib * S_i_half;

  // Compute the exponential  exp(-i*Hvib_eff*dt)  
  CMATRIX expH(sz,sz);  
  exp_matrix(expH, Hvib_eff, complex<double>(0.0, -dt) );

  // Propagator
  U = S_i_half * expH * S_half;

  

}// propagate_electronic



}// namespace libelectronic
}// namespace libdyn
}// liblibra

