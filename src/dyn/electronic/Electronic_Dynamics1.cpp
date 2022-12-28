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

void propagate_electronic_qtag2(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& Hvib_old, CMATRIX& S, CMATRIX& S_old){
/**
  Solves the generalized time-dependent Schrodinger equation:

  i*hbar*S*dc/dt = Hvib*c

  using the approach described in the QTAG paper:


  U = Z(t+dt) * exp(-i/hbar *E(t+dt) *dt) * Z(t)^+  * S(t)

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
  CMATRIX Z_old(sz, sz);
  CMATRIX E(sz, sz);
  CMATRIX E_old(sz, sz);

  // Transformation to adiabatic basis
  libmeigen::solve_eigen(Hvib, S, E, Z, 0);  // Hvib(t+dt) * Z(t+dt) = S(t+dt) * Z(t+dt) * E(t+dt)
  libmeigen::solve_eigen(Hvib_old, S_old, E_old, Z_old, 0);  // Hvib(t) * Z(t) = S(t) * Z(t) * E(t)

  // Diagonal form expH
  CMATRIX expE(sz, sz);
  for(i=0;i<sz;i++){
    double argg = - 0.5*(E.get(i,i) + E_old.get(i,i)).real()*dt;
    double _cs = cos(argg);
    double _si = sin(argg);
    expE.set(i,i, complex<double>(_cs, _si)  );
  }


  // This is a supposedly more efficient form of Coeff = (Z(t+dt) * expE(t+dt/2) * Z(t).H() * S(t)) * Coeff
  Coeff = S_old * Coeff;
  Coeff = Z_old.H() * Coeff;
  Coeff = expE * Coeff;
  Coeff = Z * Coeff;

}




/*
void propagate_electronic(double dt, CMATRIX& C, nHamiltonian& ham, int rep, int method){

 propagate_electronic(dt, C, &ham, rep, method);

}


void propagate_electronic(double dt, CMATRIX& C, nHamiltonian* ham, int rep, int method){

  if(rep==0){  // diabatic
    CMATRIX Hvib(ham->ndia, ham->ndia);  Hvib = ham->get_hvib_dia();
    CMATRIX Sdia(ham->ndia, ham->ndia);  Sdia = ham->get_ovlp_dia();

    if(method==0 || method==100){
      propagate_electronic(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }
    else if(method==1 || method==101){
      propagate_electronic_qtag(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }

  }

  else if(rep==1){  // adiabatic
    CMATRIX Hvib(ham->nadi, ham->nadi);  Hvib = ham->get_hvib_adi(); 

    if(metho==0 || method==100){
      propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs
    }

  }

}
*/

void propagate_electronic(double dt, CMATRIX& C, nHamiltonian& ham, nHamiltonian& ham_prev, int rep, int method){

 propagate_electronic(dt, C, &ham, &ham_prev, rep, method);

}

void propagate_electronic(double dt, CMATRIX& C, nHamiltonian* ham, nHamiltonian* ham_prev, int rep, int method){

  if(rep==0){  // diabatic
    CMATRIX Hvib(ham->ndia, ham->ndia);
    CMATRIX Sdia(ham->ndia, ham->ndia);

    if(method==0 || method==100){
      // Based on Lowdin transformations, using mid-point Hvib
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      propagate_electronic(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }
    else if(method==1 || method==101){
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      propagate_electronic_qtag(dt, C, Hvib, Sdia); // in this case C - diabatic coeffs
    }
    else if(method==2 || method==102){
      Hvib = ham->get_ham_dia();
      Sdia = ham->get_ovlp_dia();
      CMATRIX Hvib_old(ham->ndia, ham->ndia);   Hvib_old = ham_prev->get_ham_dia();
      CMATRIX Sdia_old(ham->ndia, ham->ndia);   Sdia_old = ham_prev->get_ovlp_dia();

      propagate_electronic_qtag2(dt, C, Hvib, Hvib_old, Sdia, Sdia_old);
    }
    else if(method==3 || method==103){
      // Using exp(S^-1 * Hvib_dia * dt)
      Hvib = 0.5 * (ham->get_hvib_dia() + ham_prev->get_hvib_dia());
      Sdia = ham->get_ovlp_dia();
      CMATRIX invS(ham->ndia, ham->ndia); 
      FullPivLU_inverse(Sdia, invS);
      Hvib = invS * Hvib;
      propagate_electronic_nonHermitian(dt, C, Hvib);
    }

  }// rep == 0 // diabatic

  else if(rep==1){  // adiabatic

    if(method==0 || method==100){
      CMATRIX Hvib(ham->nadi, ham->nadi);  Hvib = 0.5*(ham->get_hvib_adi() + ham_prev->get_hvib_adi());

      cout<<"Integration with Hvib = "; Hvib.show_matrix();
      propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs
    }
    else if(method==1 || method==101){
      // Local diabatization
      CMATRIX U_old(ham->nadi, ham->nadi);   U_old = ham_prev->get_basis_transform();
      CMATRIX U(ham->nadi, ham->nadi);   U = ham->get_basis_transform();

      CMATRIX st(ham->nadi, ham->nadi); st = U_old * U.H();
      CMATRIX st2(ham->nadi, ham->nadi); st2 = st.H() * st; 
      CMATRIX st2_half(ham->nadi, ham->nadi);
      CMATRIX st2_i_half(ham->nadi, ham->nadi);
      CMATRIX T(ham->nadi, ham->nadi);

      sqrt_matrix(st2, st2_half, st2_i_half);
      T = st * st2_i_half;

      CMATRIX Z(ham->nadi, ham->nadi);

      Z = 0.5 * (ham_prev->get_ham_adi() + T * ham->get_ham_adi() * T.H() );

      propagate_electronic_rot(dt, C, Z); 

      C = T.H() * C;

    }// method == 1
   
  }

}


/*
void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep, int method){

//  This function propagates the coefficients from C(t-dt) to C(t)

//  dt - integration timestep
//  C - the matrix of amplitudes  nstates x ntraj
//  ham - Hamiltonian at time t (or at another point)
//  rep - representation: 0 - diabatic, 1 - adiabatic
//  method - anything above 100 (but below 200), including are the same methods as normal, but with the NBRA approximation

//           for the NBRA, only the first element of ham (ham[0]) will be used. However, we assume that
//           the C matrix still has the full dimensionality, that is it consist of ntraj columns


  if(! (method >= 100 and method <200) ){
    if(C.n_cols!=ham.size()){
      cout<<"ERROR in void propagate_electronic(double dt, CMATRIX& C, \
             vector<nHamiltonian*>& ham, vector<nHamiltonian*>& ham_prev, int rep, int method): \n";
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
    int traj1 = traj;  if(method >=100 and method <200){ traj1 = 0; }
    propagate_electronic(dt, ctmp, ham[traj1], rep, method);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }
  }

}
*/

/*
void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, int rep){

  int isNBRA = 0;
  propagate_electronic(dt, C, ham, rep, isNBRA);

}
*/

void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, vector<nHamiltonian*>& ham_prev, int rep, int method){
/**
  This function propagates the coefficients from C(t-dt) to C(t)

  dt - integration timestep 
  C - the matrix of amplitudes  nstates x ntraj
  ham - Hamiltonian at time t
  ham_prev - Hamiltonian at time t - dt (or any other interval ends)
  rep - representation: 0 - diabatic, 1 - adiabatic   
  method - anything above 100 (but below 200), including are the same methods as normal, but with the NBRA approximation

           for the NBRA, only the first element of ham (ham[0]) will be used. However, we assume that 
           the C matrix still has the full dimensionality, that is it consist of ntraj columns
*/


  if(! (method >= 100 and method <200) ){
    if(C.n_cols!=ham.size()){
      cout<<"ERROR in void propagate_electronic(double dt, CMATRIX& C, \
             vector<nHamiltonian*>& ham, vector<nHamiltonian*>& ham_prev, int rep, int method): \n";
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
    int traj1 = traj;  if(method >=100 and method <200){ traj1 = 0; }
    propagate_electronic(dt, ctmp, ham[traj1], ham_prev[traj1], rep, method);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }
  }

}

/*
void propagate_electronic(double dt, CMATRIX& C, vector<nHamiltonian*>& ham, vector<nHamiltonian*>& ham_prev, int rep){

  int isNBRA = 0;
  propagate_electronic(dt, C, ham, ham_prev, rep, isNBRA);

}
*/

/*
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
*/



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

}// grid_propagator


CMATRIX vectorize_density_matrix(CMATRIX* rho){

  if(rho->n_cols != rho->n_rows){
    cout<<"Error in vectorize_density_matrix: rho matrix should be a square matrix\n";
    cout<<"Current dimensions: "<<rho->n_cols<<" "<<rho->n_rows<<endl;
    cout<<"Exiting now...\n"; 
    exit(0);
  }

  int nst = rho->n_cols; 
  int sz = nst * nst;

  CMATRIX res(sz, 1);

  int cnt = 0;
  for(int i=0; i<nst; i++){
    for(int j=0; j<nst; j++){
      res.set(cnt, 0,  rho->get(i, j) );
      cnt++; 
    }
  }

  return res;
}

CMATRIX vectorize_density_matrix(CMATRIX& rho){
  return vectorize_density_matrix(&rho);
}

CMATRIX unvectorize_density_matrix(CMATRIX& rho_vec){

  int sz = rho_vec.n_rows; 
  double n_fl = sqrt(sz);
  int n = int(n_fl);
  if( n*n != sz ){  
  
    cout<<"Error in unvectorize_density_matrix: The number of elements in the input vector should be a full square of an integer\n";
    cout<<"The current size is "<<sz<<endl;
    cout<<"Exiting...\n";
    exit(0);
  }
  
  CMATRIX res(n, n);

  int cnt = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      res.set(i, j, rho_vec.get(cnt, 0));
      cnt++; 
    }
  }
  
  return res;

}

CMATRIX make_Liouvillian(CMATRIX& ham){

  if(ham.n_cols != ham.n_rows){
    cout<<"Error in make_Liouvillian: Hamiltonian matrix should be a square matrix\n";
    cout<<"Current dimensions: "<<ham.n_cols<<" "<<ham.n_rows<<endl;
    cout<<"Exiting now...\n";
    exit(0);
  }

  int nst = ham.n_cols;
  int sz = nst * nst;

  complex<double> Lijab; 
  CMATRIX L(sz, sz);
  int ij, ab;

  ij = 0;
  for(int i=0; i<nst; i++){
    for(int j=0; j<nst; j++){
      
      ab = 0;
      for(int a=0; a<nst; a++){
        for(int b=0; b<nst; b++){
          ab++;
          
          Lijab = 0.0; 
          if(j==b){  Lijab += ham.get(i, a); }
          if(i==a){  Lijab -= std::conj(ham.get(j, b)); }

          L.set(ij, ab,  Lijab );

        }// for b
      }// for a

      ij++;

    }// for j
  }// for i

  return L;
}


}// namespace libelectronic
}// namespace libdyn
}// liblibra

