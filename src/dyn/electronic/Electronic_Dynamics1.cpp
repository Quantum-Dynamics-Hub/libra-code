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

#include "Electronic.h"

/****************************************************************************
  This file contains following functions:

  void propagate_electronic(double dt,Electronic& el,Hamiltonian& pot)

****************************************************************************/

#include "../../mmath/libmmath.h"
#include <cmath>


using namespace libmmath;
using namespace libmmath::libmeigen;


namespace libdyn{
namespace libelectronic{


using libmmath::liboperators::rotate;


void Electronic::propagate_electronic(double dt,Hamiltonian* ham){
  libdyn::libelectronic::propagate_electronic(dt,this,ham);
}

void Electronic::propagate_electronic(double dt,Hamiltonian& ham){
  libdyn::libelectronic::propagate_electronic(dt,this,&ham);
}


void propagate_electronic(double dt,Electronic* el,Hamiltonian* ham){
// This version is based on the Hamiltonian formulation of TD-SE
// This propagator is good for general Hamiltonian - diabatic of adiabatic
//
// iL = iL_qp + iL_qq + iL_pp
//
// iL_qp = sum_i,j {}
//
// ...


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


void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib){
// This version is based on the Hamiltonian formulation of TD-SE
// This propagator is good for general Hamiltonian - diabatic of adiabatic
//
// iL = iL_qp + iL_qq + iL_pp
//
// iL_qp = sum_i,j {}
//
// ...


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


void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, MATRIX& S){
//
// i*hbar*S*dc/dt = Hvib*c
//
//!!! this is good only for time-independent S matrix

  int i,j;

// Let us first diagonalize the overlap matrix S

  int sz = S.num_of_cols;
  
  MATRIX* I; I = new MATRIX(sz, sz);  I->Init_Unit_Matrix(1.0);
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);

  // Transformation to adiabatic basis
  libmmath::libmeigen::solve_eigen(sz, S, *I, *Seig, *C);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


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



void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S){
//
// i*hbar*S*dc/dt = Hvib*c 
//
//!!! this is good only for time-independent S matrix

  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Transformation to adiabatic basis
  libmmath::libmeigen::solve_eigen(sz, S, *I, *Seig, *C);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


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
  libmmath::libmeigen::solve_eigen(sz, *Hvib_eff, *I, *Seig, *C);  // Hvib_eff * C = I * C * Seig  ==>  Hvib = C * Seig * C.H()


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
  delete coeff;
  delete Hvib_eff;
  delete S_i_half;
  delete S_half;
  delete Seig;
  delete C;
  delete I;



}// propagate_electronic



void propagate_electronic(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S, CMATRIX& Sdot){
//
// i*hbar*S*dc/dt = Hvib*c 
//
// this is good also for time-dependent S matrix

  int i,j;
 
  // Let us first diagonalize the overlap matrix S
  int sz = S.n_cols;  
  CMATRIX* I; I = new CMATRIX(sz, sz);  I->load_identity(); //->Init_Unit_Matrix(complex<double>(1.0,0.0));
  CMATRIX* C; C = new CMATRIX(sz, sz);  *C = complex<double>(0.0, 0.0);
  CMATRIX* Seig; Seig = new CMATRIX(sz, sz);  *Seig = complex<double>(0.0,0.0);


  // Transformation to adiabatic basis
  libmmath::libmeigen::solve_eigen(sz, S, *I, *Seig, *C);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()


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
  *Hvib_eff = (*S_i_half) * (Hvib) * (*S_i_half);  // Hermitian part
  *Hvib_eff = *Hvib_eff + complex<double>(0.0, 0.5)*(*S_i_half)*(*S_i_half)*Sdot; // non-Hermitian part: (i*hbar/2) * S^{-1} * Sdot
  *coeff = (*S_half) * (*coeff);                 // now these are effective coefficients

  
  // Set up the object
  Electronic el_eff; el_eff = el;
  for(i=0; i<sz; i++){  el_eff.q[i] = coeff->get(i,0).real(); el_eff.p[i] = coeff->get(i,0).imag(); }


  // Compute the exponential  exp(-i*Hvib*dt)  
  libmmath::libmeigen::solve_eigen(sz, *Hvib_eff, *I, *Seig, *C);  // Hvib_eff * C = I * C * Seig  ==>  Hvib = C * Seig * C.H()


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
  delete coeff;
  delete Hvib_eff;
  delete S_i_half;
  delete S_half;
  delete Seig;
  delete C;
  delete I;



}// propagate_electronic




}// namespace libelectronic
}// namespace libdyn

