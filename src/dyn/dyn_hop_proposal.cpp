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
  \file dyn_hop_proposal.cpp
  \brief The file implements the functionality to compute the hopping probabilities and do the hops
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "electronic/libelectronic.h"
#include "Dynamics.h"
#include "dyn_control_params.h"


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


namespace bp = boost::python;




MATRIX hopping_probabilities_fssh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib){
/**
  \brief This function computes the surface hopping probabilities according to Tully's FSSH prescription. 
  The surface-hopping probabilities may be Boltzmann-corrected

  See more details in:
  (1) Tully, J. C. Molecular Dynamics with Electronic Transitions. J. Chem. Phys. 1990, 93, 1061�1071. - the original paper
  (2) Fabiano, E.; Keal, T. W.; Thiel, W. Implementation of Surface Hopping Molecular Dynamics Using Semiempirical Methods. Chem. Phys. 2008, 349, 334�347.
  Here, we generalized the formula, so it works equally well for both diabatic and adiabatic representations
  (3) Akimov, A. V. Libra: An Open-Source �Methodology Discovery� Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626�1649.

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
  (1) Wang, L.; Trivedi, D.; Prezhdo, O. V. Global Flux Surface Hopping Approach for Mixed Quantum-Classical Dynamics. J. Chem. Theory Comput. 2014, 10, 3598�3605.
  (2) Akimov, A. V. Libra: An Open-Source �Methodology Discovery� Library for Quantum and Classical Dynamics Simulations.
  J. Comput. Chem. 2016, 37, 1626�1649.


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
  denmat_dot = ( denmat *  Hvib.conj() - Hvib * denmat) * complex<double>(0.0, 1.0);



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

}// gfsh




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






vector<MATRIX> hop_proposal_probabilities(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
       nHamiltonian& ham, vector<MATRIX>& prev_ham_dia){

/**
  This function computes the hop probabilities for each trajectory to hop from any state to all states

  C - are assumed to be dynamically-consistent

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;


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

    Hvib = ham.children[traj]->get_hvib_adi();

    //Transform Hamiltonian to the dynamically-consistent form:
    Hvib = projectors[traj].H() * Hvib * projectors[traj];


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

      pop_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);
      g[traj] = compute_hopping_probabilities_lz(ham.children[traj], prms.rep_lz, p_traj, invM, prev_ham_dia[traj]);

    }

    else{
      cout<<"Error in tsh1: tsh_method can be -1, 0, 1, 2, or 3. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

  }// for traj

  return g;

}






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
} // else
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



}// namespace libdyn
}// liblibra

