/*********************************************************************************
* Copyright (C) 2024 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_hop_proposal_fssh3.cpp
  \brief The file implements the functionality to compute the hopping probabilities for the FSSH3 method
  it also includes the development ideas
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "electronic/libelectronic.h"
///#include "Dynamics.h"
//#include "../hamiltonian/libhamiltonian.h"
//#include "../io/libio.h"
#include "../opt/libopt.h"
#include "../math_meigen/mEigen.h"

///#include "dyn_control_params.h"
#include "dyn_hop_proposal.h"

/// liblibra namespace
namespace liblibra{


using namespace libmeigen;
using namespace libopt;

/// libdyn namespace
namespace libdyn{


namespace bp = boost::python;


MATRIX adjust_signs(MATRIX& J){
  MATRIX res(J);
  int sz = J.n_cols;

  for(int i=1; i<sz; i++){
    double scl = SIGN( J.get(i,0) * J.get(0,i) );
    if(scl==0){ scl = 1.0; }
    res.scale(-1, i, -scl);
  }

  return res;
}

MATRIX find_best_matrix(MATRIX& c_new, MATRIX& c_old, MATRIX& J){

/**
 sum_col is bugged - need to fix
*/

  MATRIX j(J);
  j = adjust_signs(J);
  MATRIX j1(j);
  double err1 = (c_new - j * c_old).T().sum_row(0, 2);
  //cout<<"Original sign\n";
//  (c_new - j * c_old).show_matrix();
//  j.show_matrix();


//  cout<<"Changed sign\n";
  j = MATRIX(J);
  j.scale(-1,0,-1.0);
  j = adjust_signs(j);
  MATRIX j2(j);
  double err2 = (c_new - j * c_old).T().sum_row(0, 2);
//  (c_new - j * c_old).show_matrix();
//  j.show_matrix();

  MATRIX J2(J);
  J2.scale(0, -1, -1.0);

  j = adjust_signs(J2);
  MATRIX j3(j);
  double err3 = (c_new - j * c_old).T().sum_row(0, 2);
//  cout<<"Original sign + row adjustment\n";
//  (c_new - j * c_old).show_matrix();
//  j.show_matrix();

//  cout<<"Changed sign + row adjustment\n";
  j = MATRIX(J2);
  j.scale(-1,0,-1.0);
  j = adjust_signs(j);
  MATRIX j4(j);
  double err4 = (c_new - j * c_old).T().sum_row(0, 2);
//  (c_new - j * c_old).show_matrix();
//  j.show_matrix();

//  cout<<"Errors: "<<err1<<" "<<err2<<" "<<err3<<" "<<err4<<endl;
  int choice = 1;
  double err = err1; j = j1; 
  if(err2<err) { j = j2; err = err2; choice = 2; }
  if(err3<err) { j = j3; err = err3; choice = 3; }
  if(err4<err) { j = j4; err = err4; choice = 4; }
//  cout<<"Choice: "<<choice<<endl;
  

  return j;
}

/*
MATRIX find_best_matrix2(MATRIX& c_new, MATRIX& c_old, MATRIX& J){
  int a, b; 
  a = b = 0; 

  //j = MATRIX(J);
//  double err = (c_new - j * c_old).T().sum_row(0, 2);

//  for()

}

*/



vector<double> hopping_probabilities_fssh3(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& denmat_old, 
               int act_state_indx, vector<double>& errors){
/**
  \brief This function computes the surface hopping probabilities according to new experimental idea

  See more details in: TBD

  \param[in] key parameters needed for this type of calculations
    - dt - integration timestep [a.u.]
    - Temperature - temperature [ K ]
    - use_boltz_factor - whether to scale the computed probabilities by a Boltzmann factor
  \param[in] denmat - [nstates x nstates] - current density matrix
  \param[in] denmat_old - [nstates x nstates] - previous density matrix
  \param[in] act_state_indx - index of the initial state

  Returns: A nstates-vector of hopping probabilities to all states from the current active state
*/

  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k, cnt;
  double sum,g_ij,argg, nrm;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;

  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);


  //====== Testing options =========

  int size_option = prms.fssh3_size_option; // 0 - N elements, only populations; 1 - N^2 elements - also coherences
  int approach_option = prms.fssh3_approach_option; // 0 - master equation (J contains hopping probabilities); 1 - kinetic approach (J contains fluxes)
  int decomp_option = prms.fssh3_decomp_option; // 0 - bdcSvd; 1 - fullPivLu; 2 - fullPivHouseholderQr; 3 - completeOrthogonalDecomposition

  int do_adjustment =  0;
  
  //=========================== Step 1: definition of vectorized densities ==================================

  int nst2;
  if(size_option==0){ nst2 = nstates; }
  else  if(size_option==1){ nst2 = nstates * nstates; }

  MATRIX rho_old(nst2, 2);
  MATRIX rho(nst2, 2);
  MATRIX drhodt(nst2, 2);

  MATRIX PP(nstates, nstates);
  MATRIX PP_old(nstates, nstates);
  MATRIX J2(nstates, nstates);
  MATRIX P_old(nstates, 1);
  MATRIX P_new(nstates, 1);
  MATRIX dPdt(nstates, 1);


  for(i=0;i<nstates;i++){
    for(j=0;j<nstates;j++){
      complex<double> rho_ij = denmat.get(i,j);
      double val = (rho_ij * std::conj(rho_ij)).real();
      PP.set(i,j, val);      
      if(j==i){ P_new.set(i, 0, rho_ij.real()); }

      rho_ij = denmat_old.get(i,j);
      val = (rho_ij * std::conj(rho_ij)).real();
      PP_old.set(i,j, val);
      if(j==i){ P_old.set(i, 0, rho_ij.real()); }

    }// for j
  }// for i

  // Use only populations
  if(size_option==0){
    for(i=0;i<nstates;i++){
      rho_old.set(i, 0, denmat_old.get(i,i).real());
      rho.set(i, 0, denmat.get(i,i).real());
      rho_old.set(i, 1, 1.0);
      rho.set(i, 1, 1.0);
    }// for i

  }// size_option = 0 

  else if(size_option==1){  // Use the full density matrix - also coherences
    for(i=0;i<nstates;i++){
      rho_old.set(i, 0, denmat_old.get(i,i).real());
      rho.set(i, 0, denmat.get(i,i).real());
      rho_old.set(i, 1, 1.0);
      rho.set(i, 1, 1.0);
    }// for i

    int shift = int(nstates * (nstates-1) /2); // how many coherences
    cnt = nstates;
    for(i=0;i<nstates;i++){
      for(j=i+1;j<nstates;j++){
        rho_old.set(cnt, 0, denmat_old.get(i,j).real());
        rho_old.set(cnt+shift, 0, denmat_old.get(i,j).imag());
        rho.set(cnt, 0, denmat.get(i,j).real());
        rho.set(cnt+shift, 0, denmat.get(i,j).imag());
        cnt++;
      }// for j>i
    }// for i

  }// size_option = 1

  //=========================== Step 2: least squares solving ==================================

  // Derivative of the density matrix
  drhodt = (1.0/dt)*(rho - rho_old);
  dPdt = (1.0/dt)*(P_new - P_old);

  // LHS matrix = A
  MATRIX A(nst2, nst2);

  // RHS matrix = b
  MATRIX b(nst2, nst2);

  if(approach_option==0){  
    A = rho_old * rho_old.T();
    b = rho_old * rho.T();  
  }
  else if(approach_option==1){  
//    A = rho_old * rho_old.T();
//    b = rho_old * drhodt.T(); 
  }
  else if(approach_option==2){
    A = rho_old * rho_old.T();
    b = drhodt * rho_old.T();
/*
    MATRIX x(rho); x = 0.5*(rho_old + rho);
    A = x * x.T();
    b = x * drhodt.T();
*/
    
  }

  // Least squares solution:
  MATRIX J(nst2, nst2);
  //vector<double> err(5,0.0);

  if(approach_option==0){
    J = run_opt(P_old, P_new, prms.fssh3_dt, prms.fssh3_max_steps, prms.fssh3_err_tol, decomp_option, errors, dt, approach_option);
    J =  transform_to_hyperbolic(J, 10.0);
    //normalize_transition_matrix(J, 1);  SHOULD NORMALIZE!
  }
  else if(approach_option==1 || approach_option==2){
    J = run_opt(P_old, P_new, prms.fssh3_dt, prms.fssh3_max_steps, prms.fssh3_err_tol, decomp_option, errors, dt, approach_option);

    int n = P_old.n_rows;
    MATRIX num_coeff(n,n);
    MATRIX denom_coeff(n,n);

    for(j=0;j<n;j++){
      double prob = 0.0;
      if(P_old.get(j,0)>0.0){ prob = -dt*dPdt.get(j,0)/P_old.get(j,0); }
      else{ prob = 0.0; }
      if(prob<0.0){ prob = 0.0; }
      if(prob>1.0){ prob = 1.0; }
     
      for(i=0;i<n;i++){
        if(j==i){  num_coeff.set(i, j, 0.0);  denom_coeff.set(i,j, 0.0); }
        else{  num_coeff.set(i, j, prob);  denom_coeff.set(i,j, 1.0);}
      } 

    }// for i

//    cout<<"P_old = \n"; P_old.show_matrix();
//    cout<<"P_new = \n"; P_new.show_matrix();
//    cout<<"Optimized A = \n";
//    J.show_matrix(); 

    J = transform_to_hyperbolic(J, 10.0);

//    cout<<"Hyperbbolic J = \n";
//    J.show_matrix();
    MATRIX Jnorm(J);
    normalize_transition_matrix(Jnorm, num_coeff, denom_coeff);
//    cout<<"num_coeff = \n"; num_coeff.show_matrix();
//    cout<<"denom_coeff = \n"; denom_coeff.show_matrix();
//    cout<<"Hopping probabilities = \n";
//    Jnorm.show_matrix();




    //normalize_transition_matrix(J);
  }

/*
//  if(decomp_option==-2){
    MATRIX a(3*nstates, nstates*nstates);
    MATRIX B(3*nstates, 1);
    MATRIX x(nstates*nstates, 1);

    for(i=0;i<nstates;i++){
      for(j=0; j<nstates; j++){
        a.set(i, i*nstates + j, denmat_old.get(j,j).real() ); 
        a.set(i+nstates, i*nstates + j, 1.0);
        a.set(i+2*nstates, j*nstates + i, 1.0);
      }
      B.set(i, 0, denmat.get(i,i).real());
      B.set(i+nstates, 0, 1.0);
      B.set(i+2*nstates, 0, 1.0);
    }    
*/
/*
// Based on fluxes
    int sz = nstates*nstates-nstates;
    int shft = nstates - 1;
    MATRIX a(nstates, sz);
    MATRIX B(nstates, 1);
    MATRIX x(sz, 1);

    for(i=0;i<nstates;i++){
      cnt = 0;
      for(j=0; j<nstates; j++){        
        if(j==i){ ;; }
        else{
          a.set(i, i*shft + cnt, 0.5*(denmat_old.get(j,j).real() + denmat.get(j,j).real() )  );
          cnt++;
        }
      }// for i
      B.set(i, 0, drhodt.get(i,0));

   }// for i


    least_squares_solve(a, x, B, decomp_option);
*/
    //a = a * a.T();
    //cout<<det(a)<<endl;

/*
    int rank;
    int is_inver;
    FullPivLU_rank_invertible(A, rank, is_inver);

    if(is_inver){ ;; }
    else{ 
     cout<<"Matrix A is not invertible\n Old matrices = \n";
     rho_old.show_matrix();
     
     MATRIX tmp(rho_old); 
     nrm = 0.0;
     for(i=0;i<nstates-1; i++){  
       if( fabs(tmp.get(i,0) - tmp.get(i+1,0))<0.01 ){  
         tmp.add(i,   0, 0.0); 
         tmp.add(i+1, 0, 0.1); 
       }
       nrm += tmp.get(i,0);
     }
     nrm += tmp.get(nstates-1,0);
     // Renormalize:
     for(i=0;i<nstates; i++){  tmp.set(i, 0, tmp.get(i,0)/nrm); }
     

     A = tmp * tmp.T(); 
     FullPivLU_rank_invertible(A, rank, is_inver);
     if(is_inver==0){ cout<<"Matrix A is still not invertible at this point\n";
       A.show_matrix();
       cout<<"Matrix det = "<<det(A)<<endl;
     }
    }

    FullPivLU_inverse(A, J);  J = b * J; 

*/

/*
    MATRIX E_new(nstates, nstates);
    MATRIX E_old(nstates, nstates);
    MATRIX U_new(nstates, nstates);
    MATRIX U_old(nstates, nstates);

    solve_eigen(PP, E_new, U_new, 0);
    solve_eigen(PP_old, E_old, U_old, 0);

    J = U_new * U_old.T();

//    cout<<"P_old=\n"; P_old.show_matrix();
//    cout<<"P_new=\n"; P_new.show_matrix();    
//    cout<<"Original J=\n"; J.show_matrix();
//    cout<<"J * P_old = \n"; (J*P_old).show_matrix();
    J = find_best_matrix(P_new, P_old, J);

    double err = (P_new - J * P_old).T().sum_row(0, 2);
    if(err>0.05){  
      cout<<"Error: "<<err<<endl; 
      cout<<"P_old=\n"; P_old.show_matrix();
      cout<<"P_new=\n"; P_new.show_matrix();
      cout<<"Original J=\n"; J.show_matrix();
      cout<<"J * P_old = \n"; (J*P_old).show_matrix();
      cout<<"Active state index = "<<act_state_indx<<endl;
      cout<<"improved J =\n"; J.show_matrix();
      cout<<"J * P_old = \n"; (J*P_old).show_matrix();

    }

    MATRIX id(nstates, nstates); id.identity();
    err = (J*J.T() - id).tr();
    if(err > 0.01){    cout<<"error = "<<err<<endl; }

*/

/*  // Worked but wrong conceptually
    least_squares_solve(PP, J2, PP_old, decomp_option);
    J2 = J2.T();

    MATRIX imJ2(nstates, nstates); imJ2 = 0.0;
    CMATRIX cJ2(J2, imJ2);
    CMATRIX J2_half(nstates, nstates);
    CMATRIX J2_ihalf(nstates, nstates);

    sqrt_matrix(cJ2, J2_half, J2_ihalf); 
    J = J2_half.real();
*/


/* // for populations
    cnt = 0;
    for(i=0;i<nstates;i++){
      for(j=0; j<nstates; j++){
        J.set(i,j, x.get(cnt,0));
        cnt++;
      }
    }
*/

/*
// for fluxes:
    cnt = 0;
    for(i=0;i<nstates;i++){
      for(j=0; j<nstates; j++){

        if(i==j){ J.set(i,j, 0.0); }
        else{
          J.set(i,j, x.get(cnt,0));
          cnt++;
        }
      }
    }
*/
//    cout<<"J = \n"; J.show_matrix(); 

/*
  }
  else if(decomp_option==-1){
    int rank;
    int is_inver;
    FullPivLU_rank_invertible(A, rank, is_inver);

    if(is_inver){    FullPivLU_inverse(A, J);  J = b.T() * J;  }
    else { J.identity(); }
  }
  else{
    least_squares_solve(A, J, b, decomp_option);
    J = J.T();
  }
*/
  MATRIX adj(nst2, 1); adj = 0.0;
  MATRIX rho_old_adj(rho_old);
  if(do_adjustment==1){        
    for(i=0;i<nstates;i++){ rho_old_adj.set(i, 0, 0.0); }
    adj = J *  rho_old_adj;
  }


  //=========================== Step 3: probabilities ==========================================
  // Now convert them the fluxes that go from node j to all other nodes to probabilities:
  // if the flux is positive - the hop happens to the state j itself, so we don't consider them
  // if the flux is negative, we count it
  i = act_state_indx;

  double P_out = 0.0;  // total probability of leaving state i
  double Pii = rho_old.get(i,0); //0.5*(rho_old.get(i,0) + rho.get(i,0));
//  cout<<"Pii = "<<Pii<<endl;
  if(Pii > 0.0){
    P_out = -(dt*drhodt.get(i,0) - adj.get(i,0))/Pii;
    if(P_out<0.0){ P_out = 0.0; }
    if(P_out>1.0){ P_out = 1.0; }
  }else { P_out = 0.0; }

//  cout<<"P_out = "<<P_out<<endl;
  
  nrm = 0.0;
  for(j=0; j<nstates;j++){
    if(approach_option==-1){ //GFSH
      if(j!=i){
        g[j] = drhodt.get(j,0);
        if(g[j] > 0.0) { nrm += g[j]; }
        else{ g[j] = 0.0; }
      }
    }  
    else if(approach_option==0){   // Option for transition probabilities    
      g[j] = J.get(j,i);
      if(g[j] > 0.0) { nrm += g[j]; }
      else{ g[j] = 0.0; }
    }
    else if(approach_option==1 || approach_option==2){ // Option for fluxes
      g[j] = J.get(j,i);
//      cout<<"g["<<j<<"] = "<<g[j]<<" ";
      if(j!=i){
        if(g[j] > 0) { nrm += g[j]; }
        else{ g[j] = 0.0; }
      }
//      cout<<"g_corr["<<j<<"] = "<<g[j]<<" ";
    }
//    cout<<endl;
  }// for j - all states

//  if(fabs(nrm)<1e-10) { nrm = 1.0; }

  // Normalize  
  for(j=0; j<nstates;j++){
    if(approach_option==0){  
      // For the populations approach, we already have probabilities
      // but let's normalize still
      g[j] = fabs(g[j]/nrm);  
    }
    else if(approach_option==1 || approach_option==2 || approach_option==-1){
      // For the fluxes approach, the staying probability is defined based on the 
      // total outflow probability and the hopping probabilities are proportional to 
      // fluxes
  
      if(j==i){  g[j] = 1.0 - P_out; }
      else{      g[j] = fabs(g[j]/nrm)*P_out;   }
      
    }
  }// for all states

//  cout<<"Final hopping probabilities\n";
//  for(i=0;i<nstates;i++){ cout<<g[i]<<" "; } cout<<endl;

  return g;


}// fssh3




}// namespace libdyn
}// liblibra
