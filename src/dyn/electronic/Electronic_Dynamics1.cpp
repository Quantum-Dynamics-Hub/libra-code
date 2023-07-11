/*********************************************************************************
* Copyright (C) 2015-2023 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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


void iL2_action(double dt, CMATRIX& Coeff, CMATRIX& Hvib, int i, int j){
// Action of iL2 on a pair of coefficients c_i and c_j

  CMATRIX exp_iL2(2, 2);
  CMATRIX C(2, 1);
  double A = dt * Hvib.get(i,j).imag();
  complex<double> one(1.0, 0.0);
  double cs = cos(A);
  double si = sin(A);

  exp_iL2.set(0, 0,  cs*one);    exp_iL2.set(0, 1,  si*one);
  exp_iL2.set(1, 0, -si*one);    exp_iL2.set(1, 1,  cs*one);

  C.set(0,0, Coeff.get(i,0));
  C.set(1,0, Coeff.get(j,0));

  C = exp_iL2 * C;

  Coeff.set(i,0, C.get(0,0));
  Coeff.set(j,0, C.get(1,0));

}

void iL3_action(double dt, CMATRIX& Coeff, CMATRIX& Hvib, int i, int j){
// Action of iL3 on a pair of coefficients c_i and c_j

  CMATRIX exp_iL3(2, 2);
  CMATRIX C(2, 1);
  double B = dt * Hvib.get(i,j).real();
  complex<double> one(1.0, 0.0);
  complex<double> eye(0.0, 1.0);
  double cs = cos(B);
  double si = sin(B);

  exp_iL3.set(0, 0,  cs*one);    exp_iL3.set(0, 1, -si*eye);
  exp_iL3.set(1, 0, -si*eye);    exp_iL3.set(1, 1,  cs*one);

  C.set(0,0, Coeff.get(i,0));
  C.set(1,0, Coeff.get(j,0));

  C = exp_iL3 * C;

  Coeff.set(i,0, C.get(0,0));
  Coeff.set(j,0, C.get(1,0));

}

void propagate_electronic_rot(double dt, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  \brief Propagate electronic DOF using sequential rotations in the amplitudes variables, not MMTS!

  Solves: i * hbar * dC/dt = Hvib * C, where C = Coeff 

  ** phases ** - it is norm-conserving
  iL^(1) = \sum_i { -i/hbar * H_ii * d/dC_i} 
 
   action:  C_i(t) ->  C(t+dt) = exp(- i/hbar * H_ii * dt ) * C(t)

  ** nonadiabatic population transfer ** - it is norm-conserving
  iL^2(2) = \sum_{i>i} {  Im(H_ij)/hbar * ( C_j d/dC_i - C_i d/dC_j )}  

   action: 
           (C_i(t) )        (C_i(t+dt))        (   cos(A)   sin(A) )  ( C_i(t) )
           (       )   ->   (         )    =   (                   )  (        )
           (C_j(t) )        (C_j(t+dt))        (  -sin(A)   cos(A) )  ( C_j(t) )
  
  A = Im(H_ij)/hbar

  ** diabatic population transfer ** - it is norm-conserving
  iL^2(3) = \sum_{i>i} {  Re(H_ij)/hbar * ( C_j d/dC_i + C_i d/dC_j )}

   action:
           (C_i(t) )        (C_i(t+dt))        (   cos(B)   -i*sin(B) )  ( C_i(t) )
           (       )   ->   (         )    =   (                      )  (        )
           (C_j(t) )        (C_j(t+dt))        (  -i*sin(B)   cos(B)  )  ( C_j(t) )
         
  B = Re(H_ij)/hbar

  \param[in] dt The integration time step (also the duration of propagation)
  \param[in,out] Coeff(nst, 1) the adiabatic amplitudes
  \param[in] Hvib(nst, nst) the vibronic Hamiltonian, must be Hermitian

  This version can have non-zero real components of the off-diagonal elements of Hvib
  as long as Hvib is kept Hermitian

  This is the Python-friendly function
*/


  int i,j;
  double ci, cj;

  double dt_half = 0.5*dt;
  complex<double> tau(0.0, -0.5*dt);  // this is actually -i*dt/2

  if(Coeff.n_cols!=1){
    cout<<"ERROR in propagate_electronic_rot_new: The number of columns in the input amplitudes vector\
          should be 1, but "<<Coeff.n_cols<<" is given \n";
    exit(0);
  }
  if(Hvib.n_cols!=Hvib.n_rows){
    cout<<"ERROR in propagate_electronic_rot_new: The number of columns and rows of the input Hamiltonian\
          should be equal to each other, but the numbers are: \n";
    cout<<"num_of_rows = "<<Hvib.n_rows<<"\n";
    cout<<"num_of_cols = "<<Hvib.n_cols<<"\n";
    exit(0);
  }
  if(Hvib.n_cols!=Coeff.n_rows){
    cout<<"ERROR in propagate_electronic_rot_new: The number of columns of the input Hamiltonian\
          should be equal to the number of rows of the amplitudes vector, but what is given: \n";
    cout<<"Coeff.num_of_rows = "<<Coeff.n_rows<<"\n";
    cout<<"Hvib.num_of_cols = "<<Hvib.n_cols<<"\n";
    exit(0);
  }


  int nstates = Coeff.n_rows;

  //------------- Phase evolution according iL^(1) ----------------
  // exp(iL^(1) * dt/2)
  for(i=0;i<nstates;i++){
      Coeff.set(i, 0,  exp(tau*Hvib.get(i,i)) * Coeff.get(i, 0) ); 
  }// for i


  //------------- Diabatic population transfer, according to iL^(3) ----------------
  // exp(iL^(3) * dt/2), forward e.g. iL_01 * iL_02 * iL_12
  for(i=0;i<nstates;i++){
    for(j=0;j<i;j++){
      iL3_action(dt_half, Coeff, Hvib, i, j);
    }// for j
  }// for i


  //------------- Non-adiabatic population transfer, according to iL^(2) ----------------
  // exp(iL^(2) * dt/2), forward e.g. iL_01 * iL_02 * iL_12
  for(i=0;i<nstates;i++){
    for(j=0;j<i;j++){
      iL2_action(dt_half, Coeff, Hvib, i, j);
    }// for j
  }// for i

  // exp(iL^(2) * dt/2), backward e.g. iL_12 * iL_02 * iL_01
  for(i=nstates-1;i>=0;i--){
    for(j=i-1;j>=0;j--){
      iL2_action(dt_half, Coeff, Hvib, i, j);
    }// for j
  }// for i


  //------------- Diabatic population transfer, according to iL^(3) ----------------
  // exp(iL^(3) * dt/2), backward e.g. iL_01 * iL_02 * iL_12
  for(i=nstates-1;i>=0;i--){
    for(j=i-1;j>=0;j--){
      iL3_action(dt_half, Coeff, Hvib, i, j);
    }// for j
  }// for i


  //------------- Phase evolution according iL^(1) ----------------
  // exp(iL^(1) * dt/2)
  for(i=0;i<nstates;i++){
      Coeff.set(i, 0,  exp(tau*Hvib.get(i,i)) * Coeff.get(i, 0) );
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
  CMATRIX I(sz, sz);  I.load_identity();
  CMATRIX C(sz, sz);  // eigenvectors
  CMATRIX Heig(sz, sz);  // eigenvalues
  CMATRIX expH(sz, sz);  

  // Compute the exponential  exp(-i*Hvib*dt)  
  libmeigen::solve_eigen(Hvib, I, Heig, C, 0);  // Hvib_eff * C = I * C * Heig  ==>  Hvib = C * Heig * C.H()

  // Diagonal form expH
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Heig.get(i,i)*dt );
    expH.set(i,i,val);
  }

  // Transform back to the original basis:
  expH = C * expH * C.H();

  // Propagation
  Coeff = expH * Coeff;
  

}// propagate_electronic_eig


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
  CMATRIX I(sz, sz);  I.load_identity();
  CMATRIX R(sz, sz); // right eigenvectors:  A * R = R * D
  CMATRIX L(sz, sz); // left eigenvectors:   L * A  = D * L <=> A.H() * L.H() = L.H() * D
  CMATRIX Heig(sz, sz); // eigenvalues - same for Hvib and HvibH, but the eigenvectors are different
  CMATRIX HvibH(sz,sz); HvibH = Hvib.H(); 
  CMATRIX expH(sz, sz);  

  // Solve the right eigenvalue problem
  libmeigen::solve_eigen(Hvib, I, Heig, R, 0);  // Hvib * R = R * Heig  

  // Solve the left eigenvalue problem
  libmeigen::solve_eigen(HvibH, I, Heig, L, 0);  // Hvib.H() * L.H() = L.H() * Heig  
  L = L.H();


  // Diagonal form expH
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Heig.get(i,i)*dt );
    expH.set(i,i,val);
  }

  // Transform back to the original basis:
  expH = R * expH * L;

  // Propagation
  Coeff = expH * Coeff;


}// propagate_electronic_nonHermitian


void propagate_electronic_eig(double dt, CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S){
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
  CMATRIX I(sz, sz);  I.load_identity();
  CMATRIX C(sz, sz);    // S eigenvectors
  CMATRIX Seig(sz, sz); // S eigenvalues

  // Transformation to adiabatic basis
  libmeigen::solve_eigen(S, I, Seig, C, 0);  // S * C = I * C * Seig  ==>  S = C * Seig * C.H()

  // Diagonal form of S^{-1/2} and S^{1/2} matrices
  CMATRIX S_i_half(sz, sz);  // S^{-1/2}
  CMATRIX S_half(sz, sz);    // S^{1/2}

  for(i=0;i<sz;i++){
    complex<double> val = std::sqrt(Seig.get(i,i));
    S_i_half.set(i,i, 1.0/val);
    S_half.set(i,i, val);
  }

  // Convert to original basis
  S_i_half = C * S_i_half * C.H();
  S_half = C * S_half * C.H();

  // Transform the Hamiltonian accordingly:
  CMATRIX Hvib_eff(sz,sz);
  CMATRIX coeff(sz,1);

  Hvib_eff = S_i_half * Hvib * S_i_half;  // Hermitian part
  coeff = S_half * Coeff;                 // now these are the effective coefficients

  
  // Compute the exponential  exp(-i*Hvib*dt)
  // here we use Seig to also store the eigenvalues of Hvib_eff  
  libmeigen::solve_eigen(Hvib_eff, I, Seig, C, 0);  // Hvib_eff * C = I * C * Seig  ==>  Hvib = C * Seig * C.H()

  // Diagonal form expH
  CMATRIX expH(sz, sz);  // 
  complex<double> one(0.0, 1.0);
  for(i=0;i<sz;i++){
    complex<double> val = std::exp(-one*Seig.get(i,i)*dt );
    expH.set(i,i,val);
  }

  // Transform back to the original basis:
  expH = C * expH * C.H();

  // Propagation
  coeff = expH * coeff;

  // Transform the coefficients back to the original representation:
  coeff = S_i_half * coeff;  // convert the effective coefficients back to original representation

  for(i=0;i<Coeff.n_elts;i++){  Coeff.set(i, 0, coeff.get(i, 0)); } // so we don't allocate new memory for Coeff!!!
 
}// propagate_electronic_eig



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
         
          Lijab = 0.0; 
          if(j==b){  Lijab += ham.get(i, a); }
          if(i==a){  Lijab -= std::conj(ham.get(j, b)); }

          L.set(ij, ab,  Lijab );
          ab++;

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

