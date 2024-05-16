/*********************************************************************************
* Copyright (C) 2019-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_hop_proposal.cpp
  \brief The file implements the functionality to compute the hopping probabilities and do the hops
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


MATRIX hopping_probabilities_fssh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  \brief This function computes the surface hopping probabilities according to Tully's FSSH prescription. 
  The surface-hopping probabilities may be Boltzmann-corrected

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition - is dynamically-consistent
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix - must be dynamically-consistent
  \param[in] dt - the interval of the time step for which we compute the hopping probability
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

*/
  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;


  int nstates = Coeff.n_rows;
  MATRIX g(nstates,nstates);
  CMATRIX denmat(nstates, nstates);   
  
  denmat = Coeff * Coeff.H();

  // Now calculate the hopping probabilities
  for(i=0;i<nstates;i++){
    sum = 0.0;
    double a_ii = denmat.get(i,i).real(); 

    for(j=0;j<nstates;j++){

      if(i!=j){ 
        /**
          dc/dt = -(i/hbar) * Hvib * c
          (dc/dt)^+ = i/hbar * c^+ * Hvib^+

          rho = c * c^+

          Then
          drho/dt = i/hbar * (rho*Hvib - Hvib*rho)

          The diagonal element:
           
          drho_ii/dt = (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai}

          Then:  P(i->*) = -(drho_ii / rho_ii ) * dt  - probability of leaving state i

          Then:  P(i->a) = - (dt/rho_ii) * Re[ (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai} ] 

          = (dt/(hbar*rho_ii)) * Im[ rho_ia * Hvib_ai - Hvib_ia * rho_ai ] 

          Or:

          P(i->j) = (dt/(hbar*rho_ii)) * Im[ rho_ij * Hvib_ji - Hvib_ij * rho_ji ] 

        */

        double imHaij = ( denmat.get(i,j) * Hvib.get(j,i) - Hvib.get(i,j) * denmat.get(j,i) ).imag(); 

        if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{
          g_ij = dt*imHaij/a_ii;  // This is a general case -

          if(use_boltz_factor){

            if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
              argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
              if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
            }

          }// if use_boltz_factor


          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g.set(i,j,g_ij);
        sum = sum + g_ij;
      }
      else{ g.set(i,j,0.0); }

    }// for j

    g.set(i,i,1.0 - sum);

  }// for i

  return g;

}// fssh



vector<double> hopping_probabilities_fssh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx){
/**
  \brief This function computes the surface hopping probabilities according to Tully's FSSH prescription. 
  The surface-hopping probabilities may be Boltzmann-corrected

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061–1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334–347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

  \param[in] key parameters needed for this type of calculations
    - dt - integration timestep [a.u.]
    - Temperature - temperature [ K ]
    - use_boltz_factor - whether to scale the computed probabilities by a Boltzmann factor
  \param[in] denmat - [nstates x nstates] - density matrix 
  \param[in] Hvib - [nstates x nstates] - vibronic Hamiltonian matrix
  \param[in] act_state_indx - index of the initial state 

  Returns: A nstates-vector of hopping probabilities to all states from the current active state

*/
  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;

  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);


//  cout<<"In hopping prob fssh\n";
  // Now calculate the hopping probabilities
  i = act_state_indx;
   
  sum = 0.0;
  double a_ii = denmat.get(i,i).real(); 
//  cout<<"a_ii = "<<a_ii<<endl;

  for(j=0;j<nstates;j++){
 
    if(i!=j){ 
      /**
        dc/dt = -(i/hbar) * Hvib * c
        (dc/dt)^+ = i/hbar * c^+ * Hvib^+
 
        rho = c * c^+
 
        Then
        drho/dt = i/hbar * (rho*Hvib - Hvib*rho)
 
        The diagonal element:
         
        drho_ii/dt = (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai}
 
        Then:  P(i->*) = -(drho_ii / rho_ii ) * dt  - probability of leaving state i
 
        Then:  P(i->a) = - (dt/rho_ii) * Re[ (i/hbar) *sum_a { rho_ia * Hvib_ai - Hvib_ia * rho_ai} ] 
 
        = (dt/(hbar*rho_ii)) * Im[ rho_ia * Hvib_ai - Hvib_ia * rho_ai ] 
 
        Or:
 
        P(i->j) = (dt/(hbar*rho_ii)) * Im[ rho_ij * Hvib_ji - Hvib_ij * rho_ji ] 
 
      */
 
      double imHaij = ( denmat.get(i,j) * Hvib.get(j,i) - Hvib.get(i,j) * denmat.get(j,i) ).imag(); 
//      cout<<"rho_ij = "<<denmat.get(i,j)<<"  Hvib_ij = "<<Hvib.get(i,j)<<endl;
 
      if(a_ii<1e-8){ g_ij = 0.0; }  // avoid division by zero
      else{
        g_ij = dt*imHaij/a_ii;  // This is a general case -
 
        if(use_boltz_factor){
 
          if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
            argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
            if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
          }
 
        }// if use_boltz_factor
 
 
        if(g_ij<0.0){  g_ij = 0.0; }
 
      }// else
 
      g[j] = g_ij;
      sum = sum + g_ij;
    }
    else{ g[j] = 0.0; }
 
  }// for j
 
  g[i] = 1.0 - sum;
 
  return g;
 
}// fssh





MATRIX hopping_probabilities_gfsh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib){

/**
  \brief Compute the GFSH surface hopping probabilities for a single trajectory

  \param[in] Coeff - [ndia x 1] or a [nadi x 1] matrix of electronic basis states amplitudes in a superposition
  \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
  \param[in] dt - the interval of the time step for which we compute the hopping probability
  \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
  \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1

  Returns: A matrix with the hopping probabilities between all pairs of states is returned

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.


*/

  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;


  int nstates = Coeff.n_rows;
  MATRIX g(nstates,nstates);
  CMATRIX denmat(nstates, nstates);   

  //denmat = (Coeff * Coeff.H() ).conj();
  denmat = Coeff * Coeff.H();

  //CMATRIX* denmat_dot; denmat_dot = new CMATRIX(nstates, nstates);   
  CMATRIX denmat_dot(nstates, nstates);
  denmat_dot = ( denmat *  Hvib.H() - Hvib * denmat) * complex<double>(0.0, 1.0);



  // compute a_kk and a_dot_kk
  vector<double> a(nstates,0.0);
  vector<double> a_dot(nstates,0.0);
  double norm = 0.0; // normalization factor

  for(i=0;i<nstates;i++){
    a[i] = denmat.get(i,i).real();
    a_dot[i] = denmat_dot.get(i,i).real();

    if(a_dot[i]<0.0){ norm += a_dot[i]; } // total rate of population decrease in all decaying states

  }// for i


  // Now calculate the hopping probabilities
  for(i=0;i<nstates;i++){       
    double sumg = 0.0;

    for(j=0;j<nstates;j++){

      if(j!=i){  // off-diagonal = probabilities to hop to other states

        if(a[i]<1e-12){  g.set(i,j,0.0); }  // since the initial population is almost zero, so no need for hops
        else{

          if( fabs(norm) > 1e-12 ){
            g.set(i,j,  dt*(a_dot[j]/a[i]) * a_dot[i] / norm);  
          }
 
          if(g.get(i,j)<0.0){  // since norm is negative, than this condition means that a_dot[i] and a_dot[j] have same signs
                               // which is bad - so no transitions are assigned
            g.set(i,j,0.0);
          }
          else{  // here we have opposite signs of a_dot[i] and a_dot[j], but this is not enough yet
            if(a_dot[i]<0.0 & a_dot[j]>0.0){ ;; } // this is out transition probability, but it is already computed
            else{  g.set(i,j,0.0); } // wrong transition
          }

        }// a[i]>1e-12

          if(use_boltz_factor){

            if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
              argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
              if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
            }

          }// if use_boltz_factor


        sumg += g.get(i,j);

      }
    }// for j

    g.set(i,i, 1.0 - sumg);  // probability to stay in state i

    if(g.get(i,i)<0.0){  g.set(i,i, 0.0); }

  }// for i

//  delete denmat;
//  delete denmat_dot;
  return g;

}// gfsh



vector<double> hopping_probabilities_gfsh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx){
/**
  \brief Compute the GFSH surface hopping probabilities for a single trajectory

  Abbreviation: GFSH - global flux surface hopping
  References: 
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598–3605.
  (2) Akimov, A. V. Libra: An Open-Source “Methodology Discovery” Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626–1649.

  \param[in] key parameters needed for this type of calculations
    - dt - integration timestep [a.u.]
    - Temperature - temperature [ K ]
    - use_boltz_factor - whether to scale the computed probabilities by a Boltzmann factor
  \param[in] denmat - [nstates x nstates] - density matrix 
  \param[in] Hvib - [nstates x nstates] - vibronic Hamiltonian matrix
  \param[in] act_state_indx - index of the initial state 

  Returns: A nstates-vector of hopping probabilities to all states from the current active state

*/

  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;


  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);

  CMATRIX denmat_dot(nstates, nstates);
  denmat_dot = ( denmat *  Hvib.H() - Hvib * denmat) * complex<double>(0.0, 1.0);


  // compute a_kk and a_dot_kk
  vector<double> a(nstates,0.0);
  vector<double> a_dot(nstates,0.0);
  double norm = 0.0; // normalization factor

  for(i=0;i<nstates;i++){
    a[i] = denmat.get(i,i).real();
    a_dot[i] = denmat_dot.get(i,i).real();

    if(a_dot[i]<0.0){ norm += a_dot[i]; } // total rate of population decrease in all decaying states

  }// for i


  // Now calculate the hopping probabilities
  i = act_state_indx;
  double sumg = 0.0;

  for(j=0;j<nstates;j++){

    if(j!=i){  // off-diagonal = probabilities to hop to other states

      if(a[i]<1e-12){  g[j] = 0.0; }  // since the initial population is almost zero, so no need for hops
      else{

        if( fabs(norm) > 1e-12 ){ g[j] = dt*(a_dot[j]/a[i]) * a_dot[i] / norm;  }

        // since norm is negative, than this condition means that a_dot[i] and a_dot[j] have same signs
        // which is bad - so no transitions are assigned
        if(g[j]<0.0){  g[j] = 0.0;   }

        // here we have opposite signs of a_dot[i] and a_dot[j], but this is not enough yet
        else{  
          if(a_dot[i]<0.0 && a_dot[j]>0.0){ ;; } // this is out transition probability, but it is already computed
          else{  g[j] = 0.0; } // wrong transition
        }
      }// a[i]>1e-12

        if(use_boltz_factor){
          if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
            argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
            if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
          }
        }// if use_boltz_factor

      sumg += g[j];
    }
  }// for j

  g[i] = 1.0 - sumg;  // probability to stay in state i

  if(g[i]<0.0){  g[i] = 0.0; }

  return g;

}// gfsh




vector<double> hopping_probabilities_fssh2(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& denmat_old, int act_state_indx){
/**
  \brief This function computes the surface hopping probabilities according to Leonardo Araujo's hopping probability 
  reformulation of FSSH recipe - this one does not need nonadiabatic couplings!

  See more details in: TBD


  P_{m->n}(t,t+dt) = min{ P_{m,out}(t,t+dt), (rho_{nn}(t+dt) - rho_{nn}(t))/rho_{mm}(t) }, 

  where 

  P_{m, out}(t, t+dt) = sum_{n=0, n\neqm}^{nstates} (  (rho_{nn}(t+dt) - rho{nn}(t)) /rho_{mm}(t) )

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
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;

  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);


  // Now calculate the hopping probabilities
  i = act_state_indx;

  sum = 0.0;
  double a_ii_old = denmat_old.get(i,i).real();

  for(j=0;j<nstates;j++){

    if(i!=j){
      double a_jj_old = denmat_old.get(j,j).real();
      double a_jj     = denmat.get(j,j).real();

      if(a_ii_old < 1e-8){ g_ij = 0.0; }  // avoid division by zero
      else{
        g_ij = (a_jj - a_jj_old)/a_ii_old;  // This is a general case 

        if(g_ij<0.0){  g_ij = 0.0; }
      }// else

      g[j] = g_ij;
      sum += g_ij;
    }
    else{ g[j] = 0.0; }

  }// for j

  g[i] = 1.0 - sum;

  return g;

}// fssh2


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




MATRIX hopping_probabilities_mssh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib){
/**
   \brief Compute the MSSH surface hopping probabilities scaled by Boltzmann factor
   This is the version taking the minimal amount of input information
   \param[in] Coeff The amplitudes of different basis states in the coherent superposition. This matrix is assumed to be
   a column-vector, so of the size N x 1, where N is the number of basis excited states
   \param[in] Hvib - [ndia x ndia] or a [nadi x nadi] vibronic Hamiltonian matrix
   \param[in] use_boltz_factor A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme
   \param[in] T nuclear temperature - used in the Boltzmann factor - only if use_boltz_factor is set to 1
   The function returns the matrix g with hopping probabilities
   Abbreviation: MSSH - Markov state surface hopping
   References:
   (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.
*/

  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;

  int nstates = Coeff.n_rows;
  MATRIX g(nstates,nstates);
  CMATRIX denmat(nstates, nstates);   


  double norm; norm = (Coeff.H() * Coeff).get(0,0).real();  // <- this is the norm <PSI|PSI>

  // Calculate the hopping probabilities
  for(int j=0;j<nstates;j++){
    double gjj = (std::conj(Coeff.get(j)) * Coeff.get(j)).real()/norm; // c_j^* * c_j

    for(int i=0;i<nstates;i++){ g.set(i,j,gjj);}
  }

  if(use_boltz_factor){
    if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
      argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
      if(argg<500.0){ g_ij = g_ij * std::exp(argg); }
    }
  }// if use_boltz_factor


  return g;

}


vector<double> hopping_probabilities_mssh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx){
/**
   \brief Compute the MSSH surface hopping probabilities scaled by Boltzmann factor
   This is the version taking the minimal amount of input information

   The function returns the matrix g with hopping probabilities
   Abbreviation: MSSH - Markov state surface hopping

   References:
   (1) Akimov, A. V.; Trivedi, D.; Wang, L.; Prezhdo, O. V. Analysis of the Trajectory Surface Hopping Method from the Markov State Model Perspective. J. Phys. Soc. Jpn. 2015, 84, 094002.

   \param[in] key parameters needed for this type of calculations
     - dt - integration timestep [a.u.]
     - Temperature - temperature [ K ]
     - use_boltz_factor - whether to scale the computed probabilities by a Boltzmann factor
   \param[in] denmat - [nstates x nstates] - density matrix 
   \param[in] Hvib - [nstates x nstates] - vibronic Hamiltonian matrix
   \param[in] act_state_indx - index of the initial state 

   Returns: A nstates-vector of hopping probabilities to all states from the current active state

*/

  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k;
  double sum,g_ij,argg;

  double dt = prms.dt;
  double T = prms.Temperature;
  int use_boltz_factor = prms.use_boltz_factor;

  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);

  double norm; norm = denmat.tr().real();  // <- this is the norm <PSI|PSI>

  i = act_state_indx; 

  // Calculate the hopping probabilities
  double summ = 0.0;
  for(int j=0;j<nstates;j++){

    if(j!=i){    
      double gjj = denmat.get(j,j).real()/norm; // |c_j|^2
    
      if(use_boltz_factor){
        if(Hvib.get(j,j).real() > Hvib.get(i,i).real()){
          argg = -(Hvib.get(j,j).real() - Hvib.get(i,i).real())/(kb*T);        
          if(argg<500.0){ gjj = gjj * std::exp(argg); }
        }
      }// if use_boltz_factor
          
      g[j] = gjj;
      summ += gjj;
    }
    else{  g[j] = 0.0; }
  }// for j

  g[i] = 1.0 - summ;

  return g;

}



//MATRIX compute_hopping_probabilities_lz(nHamiltonian* ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia){
vector<double> hopping_probabilities_lz(nHamiltonian* ham, nHamiltonian* ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM){
/**
  \brief This function computes the surface hopping probabilities according to Landau-Zener formula.

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061<96>1071. - the original paper
  (2) Belyaev, A. K.; Lebedev, O. V. Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
  Phys. Rev. A 2011, 84, 014701

  \param[in] ham Is the nHamiltonian object that organizes energy-related calculations
  \param[in] ham_prev Is the nHamiltonian object of previous step or geometry - needed to locate the diabatic crossing
  \param[in] act_state_indx - index of the current adiabatic state
  \param[in] rep index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic based on the min of diabatic 
  gap, 2 - adiabatic, based on the switch of the NAC sign
  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses

  Returns:
      A nstates-vector of hopping probabilities to all states from the current active state

*/
  int i,j,k;
  int nstates;
  int ndof = ham->nnucl;

  // Determine the dimensions
  if(rep==0){ nstates = ham->ndia;  }
  else if(rep==1 || rep==2){  nstates = ham->nadi;  }

//  MATRIX g(nstates,nstates);
  vector<double> g(nstates, 0.0);

  i = act_state_indx;

  /** Diabatic representation:

   p_{i->j} = exp( - 2 * pi * H_ij^2 / (hbar*v*(|H'_ii - H'_jj| )) )
 
  */
  
  if(rep==0){

    // Diabatic Ham matrices
    MATRIX ham_dia(nstates,nstates);
    MATRIX ham_dia_prev(nstates, nstates);
    ham_dia = ham->get_ham_dia().real();
    ham_dia_prev = ham_prev->get_ham_dia().real();

    // Get denominator
    MATRIX dham(nstates,nstates); // H'_ii * v
    MATRIX dHv(nstates,nstates); // |H'_ii - H'_jj| * v

    for(k=0;k<ndof;k++){
      dham = ham->get_d1ham_dia(k).real() * p.get(k,0) * invM.get(k,0);

      for(j=0;j<nstates;j++){
        double dmod = fabs(dham.get(i,i) - dham.get(j,j));
        dHv.add(i,j, dmod);
        dHv.add(j,i, dmod);
      }// for j
    }// for k


    //========== Now calculate the hopping probabilities ===========
    
    // Off-diagonal elements - only, if we meet the gap minimum
    for(j=0; j<nstates; j++){

      if(j!=i){
        double dh = ham_dia.get(i,i) - ham_dia.get(j,j);
        double dh_prev = ham_dia_prev.get(i,i) - ham_dia_prev.get(j,j);

        if(dh * dh_prev < 0.0){

          double h_ij = ham_dia.get(i,j);
          double g_ij = exp(-2.0*M_PI*h_ij*h_ij / dHv.get(i,j) );

          g[j] = g_ij;
          //g.set(i,j,g_ij);
          //g.set(j,i,g_ij);

        }// only if minimum is located 
        else{ g[j] = 0.0; }

      }// j!=i
    }// for j
  }// rep == 0 diabatic


  // Adiabatic representation
  else if(rep==1 || rep==2){

    // Update adiabatic NACs - assume those are updated 
    //ham->compute_nac_adi(p, invM);

    // Now calculate the hopping probabilities
    MATRIX ham_dia(nstates,nstates);
    MATRIX ham_adi(nstates,nstates);
    MATRIX nac_adi(nstates,nstates);
    MATRIX ham_adi_prev(nstates,nstates);
    MATRIX ham_dia_prev(nstates,nstates);
    MATRIX nac_adi_prev(nstates,nstates);

    ham_dia = ham->get_ham_dia().real();
    ham_adi = ham->get_ham_adi().real();
    nac_adi = ham->get_nac_adi().real();
    ham_adi_prev = ham_prev->get_ham_adi().real();
    ham_dia_prev = ham_prev->get_ham_dia().real();
    nac_adi_prev = ham_prev->get_nac_adi().real();

    // Off-diagonal elements
    for(j=0;j<nstates;j++){
      if(j!=i){

        int nonzero = 0;
        if(rep==1){
          double dh = ham_dia.get(i,i) - ham_dia.get(j,j);
          double dh_prev = ham_dia_prev.get(i,i) - ham_dia_prev.get(j,j);
          if( dh * dh_prev < 0.0){ nonzero = 1; }
        }
        else if(rep==2){
          double nac = nac_adi.get(i,j);
          double nac_prev = nac_adi_prev.get(i,j);
          if( nac * nac_prev < 0.0){ nonzero = 1; }
        }

        if(nonzero==1){  // This is the condition of intersection of diabatic surfaces, alternatively, we 
                         // could be looking to locate the adiabatic surfaces minimum
          double g_ij = 0.0;
          double Z_ij = fabs(ham_adi.get(i,i) - ham_adi.get(j,j));
          double nac_ij = fabs(nac_adi.get(i,j));
          if(nac_ij > 0.0){ g_ij = exp(-0.25*M_PI*Z_ij/nac_ij);     }         
          g[j] = g_ij;
        }
        else{ g[j] = 0.0; }

      }// j!=i
    }// for j
  }// rep == 1 adiabatic

  // Diagonal elements - common for both reps
  double sum = 0.0;
  for(j=0;j<nstates;j++){
    if(j!=i){  sum += g[j];   }
  }// for i
  g[i] = 1.0 - sum;

  return g;

}// lz


//MATRIX compute_hopping_probabilities_lz(nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia){
vector<double> hopping_probabilities_lz(nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM){
/**
  \brief This function computes the surface hopping probabilities according to Landau-Zener formula.

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061<96>1071. - the original paper
  (2) Belyaev, A. K.; Lebedev, O. V. Nonadiabatic nuclear dynamics of atomic collisions based on branching classical trajectories
  Phys. Rev. A 2011, 84, 014701

  \param[in] ham Is the nHamiltonian object that organizes energy-related calculations
  \param[in] rep index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic
  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses
  \param[in] prev_ham_dia [nstates x nstates] the matrix of previous diabatic Hamiltonian - needed to locate the diabatic crossing

  Returns:
      A matrix with the hopping probabilities between all pairs of states is returned
      Convention: P(i,j) is the probability for the i->j transition

*/

  return hopping_probabilities_lz(&ham, &ham_prev, act_state_indx, rep, p, invM);

}

vector<double> hopping_probabilities_zn(nHamiltonian* ham, nHamiltonian* ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM){
/**
  \brief This function computes the surface hopping probabilities according to Zhu-Nakamura generalized to multiple dimensons 
  according to formula of Yu et al. See my Chapter, Eqs. 3.89 - 3.91

  See more details in:

  \param[in] ham Is the nHamiltonian object that organizes energy-related calculations
  \param[in] ham_prev Is the nHamiltonian object of previous step or geometry - needed to locate the diabatic crossing
  \param[in] act_state_indx - index of the current adiabatic state
  \param[in] rep index selecting the type of representation to be used: 0 - diabatic, 1 - adiabatic based on the min of diabatic
  gap, 2 - adiabatic, based on the switch of the NAC sign
  \param[in] p [ndof x 1] or [ndof x ntraj] matrix of nuclear momenta
  \param[in] invM [ndof x 1] matrix of inverse nuclear masses

  Returns:
      A nstates-vector of hopping probabilities to all states from the current active state

*/
  int i,j,k;

  int ndof = ham->nnucl;
  int nstates = ham->ndia; 
//  int ntraj = ham->children.size();

  vector<double> g(nstates, 0.0);

  i = act_state_indx;

  /** Diabatic representation:

   p_{i->j} = exp( - (pi / 4.0 * |a|) * sqrt(  2 / (b^2 + sqrt( | b^4  + sign(F_i * F_j) |)) ) )

  */

  // Diabatic Ham matrices
  MATRIX ham_dia(nstates,nstates);
  MATRIX ham_dia_prev(nstates, nstates);
  ham_dia = ham->get_ham_dia().real();
  ham_dia_prev = ham_prev->get_ham_dia().real();

  MATRIX F(ndof, nstates);
  MATRIX Fi(ndof, 1);
  MATRIX Fj(ndof, 1);
  MATRIX tmp(ndof, 1);
  vector<int> st(1, ham->id);

  F = ham->all_forces_adi(st).real().T();
  Fi = F.col(i);

  // Off-diagonal elements - only, if we meet the gap minimum
  for(j=0; j<nstates; j++){

    if(j!=i){
      double dh = ham_dia.get(i,i) - ham_dia.get(j,j);
      double dh_prev = ham_dia_prev.get(i,i) - ham_dia_prev.get(j,j);

      if(dh * dh_prev < 0.0){
        Fj = F.col(j);

        double sgn = (Fi.T() * Fj).get(0);
        sgn = SIGN(sgn);
        
        // |Fi * Fj| = sqrt( |Fi^T * iM * F_j |)
        tmp.dot_product( Fi,  Fj);
        tmp.dot_product( tmp, invM);
        double x = sqrt(fabs(tmp.sum()));

        // |F_i - F_j | = sqrt( (F_i - F_j)^T * iM * (F_i - F_j) )
        Fj = Fi - Fj;
        tmp.dot_product( Fj, Fj );
        tmp.dot_product( tmp, invM);
        double y = sqrt( fabs(tmp.sum() ) );

        double h = fabs(ham_dia.get(i,j));

        if(h > 0.0 && x > 0.0){
          double a2 = 0.0625 * x * y / (h * h * h);
          double b2 = 0.5 * y/(x * h) ;
          double g_ij = exp( -(0.25*M_PI / sqrt(a2) ) * sqrt(2.0 / ( b2 + sqrt(b2*b2 + sgn )  )) );

          g[j] = g_ij;
        }
        else{ g[j] = 0.0; }

      }// only if minimum is located
      else{ g[j] = 0.0; }

    }// j!=i
  }// for j

  // Diagonal elements - common for both reps
  double sum = 0.0;
  for(j=0;j<nstates;j++){
    if(j!=i){  sum += g[j];   }
  }// for i
  g[i] = 1.0 - sum;

  return g;

}

vector<double> hopping_probabilities_zn(nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM){

  return hopping_probabilities_zn(&ham, &ham_prev, act_state_indx, rep, p, invM);

}




vector<double> hopping_probabilities_mash(dyn_control_params& prms, CMATRIX& denmat){

  int nstates = denmat.n_rows;
  vector<double> g(nstates, 0.0);

  // In MASH, we set the hopping probability to 1.0 for the state with maximal density
  int max_dens_indx = 0;
  double max_dens = denmat.get(0,0).real();
  for(int i=1; i<nstates; i++){
    double dens = denmat.get(i,i).real();
    if(dens > max_dens){  max_dens_indx = i; max_dens = dens; }
  }

  g[max_dens_indx] = 1.0;
  
  return g;
}



/*
vector<MATRIX> hop_proposal_probabilities(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
       nHamiltonian& ham, vector<MATRIX>& prev_ham_dia){
*/
vector<MATRIX> hop_proposal_probabilities(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, 
       nHamiltonian& ham, vector<MATRIX>& prev_ham_dia){
/**
  This function computes the hop probabilities for each trajectory to hop from any state to all states

  C - are assumed to be dynamically-consistent

  TO BE DEPRECATED

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;
  int isNBRA = prms.isNBRA;


  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);


  vector<MATRIX> g(ntraj, MATRIX(nst,nst)); /// the matrices of hopping probability
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);
  CMATRIX Hvib(nst, nst);
  vector<int> fstates(ntraj,0); 

  //============== Begin the TSH part ===================  
  // Proposed hops probabilities

  for(traj=0; traj<ntraj; traj++){

    nucl_stenc_y[0] = traj;
    el_stenc_y[0] = traj;
    full_id[1] = traj;

    pop_submatrix(C, coeff, el_stenc_x, el_stenc_y);

    if(isNBRA==1){
    // Only compute the Hvib for one traj, traj==0
      if(traj==0){
      Hvib = ham.children[traj]->get_hvib_adi();

      //Transform Hamiltonian to the dynamically-consistent form:
      //Hvib = projectors[traj].H() * Hvib * projectors[traj];
      }
    }

    else{
    // Compute the Hvib for all traj
    Hvib = ham.children[traj]->get_hvib_adi();

    //Transform Hamiltonian to the dynamically-consistent form:
    //Hvib = projectors[traj].H() * Hvib * projectors[traj];
    }

    if(prms.tsh_method == 0){ // FSSH

      g[traj] = hopping_probabilities_fssh(prms, coeff, Hvib);

    }
    else if(prms.tsh_method == 1){ // GFSH

      g[traj] = hopping_probabilities_gfsh(prms, coeff, Hvib);

    }
    else if(prms.tsh_method == 2){ // MSSH

      g[traj] = hopping_probabilities_mssh(prms, coeff, Hvib);

    }
    else if(prms.tsh_method == 3){ // LZ

//      pop_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);
//      g[traj] = compute_hopping_probabilities_lz(ham.children[traj], prms.rep_lz, p_traj, invM, prev_ham_dia[traj]);

    }

    else{
      cout<<"Error in tsh1: tsh_method can be -1, 0, 1, 2, or 3. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

  }// for traj

  return g;

}


vector< vector<double> > hop_proposal_probabilities(dyn_control_params& prms, dyn_variables& dyn_var, 
nHamiltonian& ham, nHamiltonian& ham_prev){
//vector<MATRIX>& prev_ham_dia){
/**

  This function computes the hop probabilities for each trajectory to hop from the current state to all states

  Returns: array of the ntraj x nstates shape
*/

  int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;
  int nst = dyn_var.nadi;
  if(prms.rep_tdse==0 || prms.rep_tdse==2){ nst = dyn_var.ndia; }

  int isNBRA = prms.isNBRA;

  int traj, dof, i;

  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);


  vector< vector<double> > g(ntraj, vector<double>(nst,0.0) ); /// the hopping probability for all trajectories
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);
  CMATRIX Hvib(nst, nst);
  vector<int> fstates(ntraj,0); 

  //============== Begin the TSH part ===================  
  // Proposed hops probabilities
  for(traj=0; traj<ntraj; traj++){

    CMATRIX& dm = *dyn_var.dm_adi[traj];
    if(prms.rep_tdse==0 || prms.rep_tdse==2){ dm = *dyn_var.dm_dia[traj]; }

    nucl_stenc_y[0] = traj;
    el_stenc_y[0] = traj;
    full_id[1] = traj;

    //pop_submatrix(C, coeff, el_stenc_x, el_stenc_y);

    if(isNBRA==1){
      // Only compute the Hvib for one traj, traj==0
      if(traj==0){   Hvib = ham.children[traj]->get_hvib_adi();   }
    }

    else{
      // Compute the Hvib for all traj
      Hvib = ham.children[traj]->get_hvib_adi();
    }

    if(prms.tsh_method == 0){ // FSSH

      g[traj] = hopping_probabilities_fssh(prms, dm, Hvib, dyn_var.act_states[traj]);

    }
    else if(prms.tsh_method == 1){ // GFSH

      g[traj] = hopping_probabilities_gfsh(prms, dm, Hvib, dyn_var.act_states[traj]);

    }
    else if(prms.tsh_method == 2){ // MSSH

      g[traj] = hopping_probabilities_mssh(prms, dm, Hvib, dyn_var.act_states[traj]);

    }
    else if(prms.tsh_method == 3){ // LZ

      pop_submatrix(*dyn_var.p, p_traj, nucl_stenc_x, nucl_stenc_y);
      g[traj] = hopping_probabilities_lz(ham.children[traj], ham_prev.children[traj], dyn_var.act_states[traj], 
                prms.rep_lz, p_traj, *dyn_var.iM);

    }
    else if(prms.tsh_method == 4){ // ZN

      pop_submatrix(*dyn_var.p, p_traj, nucl_stenc_x, nucl_stenc_y);
      g[traj] = hopping_probabilities_zn(ham.children[traj], ham_prev.children[traj], dyn_var.act_states[traj],
                prms.rep_lz, p_traj, *dyn_var.iM);

    }
 
    else if(prms.tsh_method == 5){ // DISH, but handled differently
    }

    else if(prms.tsh_method == 6){ // MASH 
      g[traj] = hopping_probabilities_mash(prms, dm);
    }

    else if(prms.tsh_method == 7){ // FSSH2

      CMATRIX& dm_prev = *dyn_var.dm_adi_prev[traj];
      if(prms.rep_tdse==0 || prms.rep_tdse==2){ dm = *dyn_var.dm_dia_prev[traj]; }

      g[traj] = hopping_probabilities_fssh2(prms, dm, dm_prev, dyn_var.act_states[traj]);
    }
    else if(prms.tsh_method == 8){ // FSSH3

      CMATRIX& dm_prev = *dyn_var.dm_adi_prev[traj];
      if(prms.rep_tdse==0 || prms.rep_tdse==2){ dm = *dyn_var.dm_dia_prev[traj]; }

      //cout<<"Coordinates\n";    dyn_var.q->show_matrix();

      CMATRIX& U = *dyn_var.proj_adi[traj];
      CMATRIX dm_prev_trans(*dyn_var.dm_adi_prev[traj]);
  
      // |psi'> = |psi> T =>  Hvib' = <psi'|H\hat|psi'> = T_new.H() * Hvib * T_new;
      // |Psi> = |psi> C = |psi'> C' = |psi> T C', so C = T C', C' = T^+ C => rho' = C' C'^+ = T^+ C C^+ T = T^+ rho T
      //  dm_adi_prev is the DM computed before the current step basis update, so it would transform as
      //  dm_adi_prev -> T^+ dm_adi_prev * T
      //dm_prev_trans = U.H() * dm_prev_trans * U;// transformation of the previousely-computed density matrix according to 
                                                // new ordering found at the current time-step (which is already reflected in the current DM).
      // However, for the GFSH option to work properly, we should not be using the transformation of the DM here      
      g[traj] = hopping_probabilities_fssh3(prms, dm, dm_prev_trans, dyn_var.act_states[traj], dyn_var.fssh3_errors[traj]);
    }

    else{
      cout<<"Error in tsh1: tsh_method can be -1, 0, 1, 2, 3, 4, 5, 6, 7, or 8. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

  }// for traj

  return g;

}






int hop(vector<double>& prob, double ksi){
/** 
  \brief Attempts a stochastic hop: basically select one of the intervals in proportion
   to their length
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/
  int i;
  int nstates = prob.size();
  double left, right; left = right = 0.0;
  int finstate = -1;
  

  // To avoid problems, lets renormalize the hopping probabilities
  double nrm = 0.0;
  if(ksi >= 1){ finstate = nstates; } // hops to highest state if ksi = 1 to avoid rounding errors
  else{
    for(i=0;i<nstates;i++){  nrm += prob[i];  }

    if(nrm>0.0){
      for(i=0;i<nstates;i++){    
        if(i==0){left = 0.0; right = prob[i]/nrm; }
        else{  left = right; right = right + prob[i]/nrm; }

        if((left<=ksi) && (ksi<=right)){  finstate = i;  }    
      } // for

    } // if nrm
    else{  finstate = initstate; }  // probability to hop to any other states is zero
                                    // so stay on the original state
                                  
    if(finstate==-1){
      std::cout<<"Something is wrong in the hop(...) function\nExiting now...\n";
      exit(0);
    }
  } // else  ksi < 1.0

  return finstate;

}// hop





int hop(int initstate, MATRIX& g, double ksi){
/** 
  \brief Attempts a stochastic hop from the initial state "initstate"
  \param[in] initstate The index of the state from which we try to hop out 
  \param[in] g The hopping probabilities matrix (type MATRIX)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/
  int i;
  int nstates = g.n_cols;
  double left, right; left = right = 0.0;
  int finstate = -1;
  

  // To avoid problems, lets renormalize the hopping probabilities
  double nrm = 0.0;
  if(ksi >= 1){ finstate = nstates; } // hops to highest state if ksi = 1 to avoid rounding errors
  else{
    for(i=0;i<nstates;i++){  nrm += g.get(initstate, i);  }

    if(nrm>0.0){
      for(i=0;i<nstates;i++){    
        if(i==0){left = 0.0; right = g.get(initstate,i)/nrm; }
        else{  left = right; right = right + g.get(initstate,i)/nrm; }

        if((left<=ksi) && (ksi<=right)){  finstate = i;  }    
      } // for

    } // if nrm
    else{  finstate = initstate; }  // probability to hop to any other states is zero
                                    // so stay on the original state
                                  
    if(finstate==-1){
      std::cout<<"Something is wrong in the hop(...) function\nExiting now...\n";
      exit(0);
    }
  } // else  ksi < 1.0

  return finstate;

}// hop




int hop(int initstate, vector<double>& g, double ksi){
/** 
  \brief Attempts a stochastic hop from the initial state "initstate"
  \param[in] initstate The index of the state from which we try to hop out 
  \param[in] g The hopping probabilities from the current active state to all other states (type vector<double>)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/
  int i;
  int nstates = g.size();
  double left, right; left = right = 0.0;
  int finstate = -1;
  

  // To avoid problems, lets renormalize the hopping probabilities
  double nrm = 0.0;
  if(ksi >= 1){ finstate = nstates; } // hops to highest state if ksi = 1 to avoid rounding errors
  else{
    for(i=0;i<nstates;i++){  nrm += g[i];  }

    if(nrm>0.0){
      for(i=0;i<nstates;i++){    
        if(i==0){left = 0.0; right = g[i]/nrm; }
        else{  left = right; right = right + g[i]/nrm; }

        if((left<=ksi) && (ksi<=right)){  finstate = i;  }    
      } // for

    } // if nrm
    else{  finstate = initstate; }  // probability to hop to any other states is zero
                                    // so stay on the original state
                                  
    if(finstate==-1){
      std::cout<<"Something is wrong in the hop(...) function\nExiting now...\n";
      exit(0);
    }
  } // else  ksi < 1.0

  return finstate;

}// hop




vector<int> propose_hops(vector<MATRIX>& g, vector<int>& act_states, Random& rnd){

  //============== Compute the proposed hops =======================
  int ntraj = act_states.size();
  vector<int> fstates(ntraj,0); 

  for(int traj=0; traj<ntraj; traj++){

    double ksi = rnd.uniform(0.0,1.0);             /// generate random number 
    fstates[traj] = hop(act_states[traj], g[traj], ksi); /// Proposed hop

  }

  return fstates;

}


vector<int> propose_hops(vector< vector<double> >& g, vector<int>& act_states, Random& rnd){

  //============== Compute the proposed hops =======================
  int ntraj = act_states.size();
  vector<int> fstates(ntraj,0); 

  for(int traj=0; traj<ntraj; traj++){

    double ksi = rnd.uniform(0.0,1.0);             /// generate random number 
    fstates[traj] = hop(act_states[traj], g[traj], ksi); /// Proposed hop

  }

  return fstates;
}



}// namespace libdyn
}// liblibra

