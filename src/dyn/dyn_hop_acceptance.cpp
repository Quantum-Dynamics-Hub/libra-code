/*********************************************************************************
* Copyright (C) 2019-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_hop_acceptance.cpp
  \brief The file implements the functions to determine if to accept hops
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "electronic/libelectronic.h"
///#include "Dynamics.h"
#include "dyn_control_params.h"
#include "dyn_hop_acceptance.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;

/// libdyn namespace
namespace libdyn{


namespace bp = boost::python;


int can_rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, vector<int>& which_dofs){
/**
  Check if we can rescale momenta p along the vector t such that the total energy is conserved

*/

  int ndof = which_dofs.size();

  // According to Fabiano
  double a_ij = 0.0;
  double b_ij = 0.0;
  double c_ij = E_old - E_new;

  for(int idof = 0; idof < ndof; idof++){  

    int dof = which_dofs[idof];

    a_ij += t.get(dof, 0) * t.get(dof, 0) * invM.get(dof, 0); 
    b_ij += p.get(dof, 0) * t.get(dof, 0) * invM.get(dof, 0); 

  }// for dof
 
  a_ij *= 0.5;

  double det = b_ij*b_ij + 4.0*a_ij*c_ij;

  int res = 0;
  if(det >= 0.0){ res = 1; }

  return res;
}

int can_rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof, 0);
  for(int i = 0; i < ndof; i++){ which_dofs[i] = i; }

  return can_rescale_along_vector(E_old, E_new, p, invM, t, which_dofs);

}



void rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse, vector<int>& which_dofs){
/**
  Rescale momenta p along the vector t such that the total energy is conserved

  or reverse the momentum along the vector t

*/
  int idof, dof; 
  int ndof = which_dofs.size();

  // According to Fabiano
  double a_ij = 0.0;
  double b_ij = 0.0;
  double c_ij = E_old - E_new;
    
  for(idof=0; idof < ndof; idof++){  

    dof = which_dofs[idof];

    a_ij += t.get(dof, 0) * t.get(dof, 0) * invM.get(dof, 0); 
    b_ij += p.get(dof, 0) * t.get(dof, 0) * invM.get(dof, 0); 

  }// for dof
 
  a_ij *= 0.5;

  double det = b_ij*b_ij + 4.0*a_ij*c_ij;

  // Calculate the scaling factor and new state
  double gamma_ij = 0.0;

  if(det<0.0){  // Frustrated hops!

    if(fabs(a_ij)>1e-100){  // only consider reversals, if the couplings are sizable
      if(do_reverse){     gamma_ij = b_ij / a_ij;}
      else{ gamma_ij = 0.0;  }
    }
    else{ gamma_ij = 0.0;  }

  }
  else{    // Accepted hops!

    if(fabs(a_ij)>1e-100){  // only consider reversals, if the couplings are sizable
      if(b_ij<0){ gamma_ij = 0.5*(b_ij + sqrt(det))/a_ij; }
      else{       gamma_ij = 0.5*(b_ij - sqrt(det))/a_ij; }
    }
    else{ gamma_ij = 0.0;  }

  }

  //Rescale velocities and do the hop
  MATRIX dp(t.n_rows, t.n_cols);
  dp = gamma_ij * t;

  /// Only shift the momenta along the active dofs (`which_dofs`)
  /// p -= dp; (only for the active dofs)
  for(idof=0; idof < ndof; idof++){  
    dof = which_dofs[idof];
    p.add(dof, 0, -dp.get(idof, 0) );
  }

}


void rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof, 0);
  for(int i = 0; i < ndof; i++){ which_dofs[i] = i; }

  rescale_along_vector(E_old, E_new, p, invM, t, do_reverse, which_dofs);

}



vector<double> Boltz_quant_prob(vector<double>& E, double T){
    /**

    Computes the quantum Boltzman probability of occupying different 
    energy levels at given temperature T:  P_i = exp(-E_i/kT) / Z
    where Z = sum_i { exp(-E_i/kT)  }

    Args:
        E ( list of doubles ): energy levels [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: prob: the probability to  find a system in a discrete state i with
            energy E_i at given temperature T, considering a number of selected energy
            levels as given by the list E
    */
    int n;
 
    double kB = boltzmann/hartree; 
    double b = 1.0/(kB * T);
    int nstates = E.size();

    double Z = 0.0;  // partition function
    vector<double> prob(nstates, 0.0);

    for(n=0; n<nstates; n++){
      prob[n] = exp(-E[n]*b);
      Z += prob[n];
    }
    
    for(n=0; n<nstates; n++){
      prob[n] = prob[n] / Z;
    }

    return prob;
}



double Boltz_cl_prob(double E, double T){
    /**

    Computes the normalized classical Boltzmann probability distribution function 

    Args: 
        E ( double ): the minimum energy level [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: The probability to have kinetic energy greater than a given threshold value at
            given temperature

    See Also:
        This is essentially a Maxwell-Boltzmann distribution in the energy scale
        Used this: http://mathworld.wolfram.com/MaxwellDistribution.html

    */
 
    int n;
 
    double kB = boltzmann/hartree; 
    double b = 1.0/(kB * T);

    double Z = 0.0;  // partition function

    double x = sqrt(E * b);
    double res = 1.0;

    double D = (2.0*b/sqrt(M_PI)) * x * exp(-x);
    if( D>1.0 || D<0.0 ){
      cout<<"D = "<<D<<endl;
      exit(0);
    }

    res = D;

    return res;

}



double Boltz_cl_prob_up(double E, double T){
    /**

    Computes the classical Boltzmann probability to have kinetic energy larger than a given 
    threshold E  at temperature T. See Eq. 7 of probabilities_theory.docx.
    The present function is related to it.

    Args: 
        E ( double ): the minimum energy level [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: The probability to have kinetic energy greater than a given threshold value at
            given temperature

    See Also:
        This is essentially a Maxwell-Boltzmann distribution in the energy scale
        Used this: http://mathworld.wolfram.com/MaxwellDistribution.html

    */

    double kB = boltzmann/hartree; 
 
    double x = sqrt(E/(kB * T));

    double D = ERF(x) - sqrt(4.0/M_PI) * x * exp(-x*x);
    if( D>1.0 || D<0.0){
      cout<<"D = "<<D<<endl;
      exit(0);
    }

    double res = 1.0 - D;

    return res;

}



double HO_prob(vector<double>& E, vector<int>& qn, double T, vector<double>& prob){

    /**
    
    Probability that the oscillators are in the given vibrational states
    Multi-oscillator generalization of Eq. 10

    Args:
        E ( list of doubles ): all the energy levels present in the system [in a.u.]
        qn ( list of ints ): quantum numbers for each oscillator
        T ( double ): temperature 

    Returns:
        tuple: (res, prob), where:

            * res ( double ): the probability that the system of N oscillators is in a given
                state, defined by quantum numbers of each oscillator
            * prob ( list of N doubles ): probability with which each of N oscillators occupies 
                given vibrational state
    
    */

    int n_freqs = E.size();

    double kB = boltzmann/hartree; 
 
    double res = 1.0;
    prob.clear();

    for(int i=0; i<n_freqs; i++){

        double xi = exp(-E[i]/(kB*T));
        prob.push_back( pow(xi, qn[i]) * (1.0 - xi) );
        res *= prob[i];
    }

    return res;
}



double HO_prob_up(vector<double>& E, vector<int>& qn, double T, vector<double>& prob){

    /**
    
    Probability that the oscillators are in the given vibrational states
    Multi-oscillator generalization of Eq. 10

    Args:
        E ( list of doubles ): all the energy levels present in the system [in a.u.]
        qn ( list of ints ): quantum numbers for each oscillator
        T ( double ): temperature 

    Returns:
        tuple: (res, prob), where:

            * res ( double ): the probability that the system of N oscillators is in a given
                state, defined by quantum numbers of each oscillator
            * prob ( list of N doubles ): probability with which each of N oscillators occupies 
                given vibrational state
    
    */

    int n_freqs = E.size();

    double kB = boltzmann/hartree; 
 
    double res = 1.0;
    prob.clear();

    for(int i=0; i<n_freqs; i++){

        double xi = exp(-E[i]/(kB*T));
        prob.push_back( xi );
        res *= prob[i];
    }

    return res;
}




double boltz_factor(double E_new, double E_old, double T, int boltz_opt){
    /**
    Compute the Boltzmann scaling factor, but only if we consider a hop up in energy

    Args: 
        E_new ( double ): the energy of the proposed (new) state [units: a.u.]
        E_old ( double ): the energy of the current (old) state [units: a.u.]
        T ( double ): temperature of the bath [units: K]
        boltz_opt ( int ): the proposed hop acceptance criterion

            * 0: all hops are accepted
            * 1: hops are accepted according to the Boltzmann ratio of the final and initial states
            * 2: hops are accepted according to the classical Maxwell-Boltzmann distribution
            * 3: hops are accepted based on the quantum probability of the final state

    Returns:
        double: boltz_f: the probability of the proposed hop acceptance

    */

    double dE = (E_new - E_old);
    double boltz_f = 1.0;
    double kB = boltzmann/hartree; 
    double argg;

    if(boltz_opt==0){
        boltz_f = 1.0;
    }
     
    else if(boltz_opt==1){
        if(dE > 0.0){
            argg = dE/(kB*T);
            if( argg > 50.0){ boltz_f = 0.0;         }
            else{             boltz_f = exp(-argg);  }
        }
    }

    else if(boltz_opt==2){
        if(dE > 0.0){
            boltz_f = Boltz_cl_prob_up(dE, T);
        }
    }

    else if(boltz_opt==3){
        if(dE > 0.0){
            vector<double> energies(2);
            energies[0] = 0.0;
            energies[1] = dE;
            boltz_f = Boltz_quant_prob(energies, T)[1];
        }
    }

    return boltz_f;

}



vector<int> accept_hops(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd, 
       vector<int>& which_trajectories){
/**
  This function returns the new state indices if the corresponding transitions can be
  accepted according to given criteria

  C - is assumed to be dynamically-consistent

  options:
  0 - accept all

  10 - based on adiabatic energy
  11 - based on diabatic energy

  20 - derivative coupling vectors
  21 - difference of state-specific forces

  31 - quantum Boltzmann
  32 - Maxwell-Boltzmann
  33 - updated quantum Boltzmann

*/

  vector<int>& which_dofs = prms.quantum_dofs;

  int ndof = q.n_rows;
  int ndof_active = which_dofs.size();
  int ntraj = q.n_cols;
  int ntraj_active = which_trajectories.size();
  int nst = C.n_rows;    
  int itraj, traj, dof, i, idof;
  int isNBRA = prms.isNBRA;

  vector<int> fstates(ntraj, -1); 
  MATRIX p_tr(ndof, 1);
  CMATRIX hvib(nst, nst);
  CMATRIX nac(nst, nst);
  CMATRIX df(nst, nst);

//  cout<<"In accept_hops with algo = "<<prms.hop_acceptance_algo<<endl;
//  cout<<"ntraj_active = "<<ntraj_active<<endl;


  if(prms.hop_acceptance_algo==0){  // Just accept all the proposed hops

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      fstates[traj] = proposed_states[traj];
    }

  }// algo = 0

  else if(prms.hop_acceptance_algo==10){  // Just based on the adiabatic energy levels

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];

      if(old_st != new_st){

        p_tr = p.col(traj);
        double T_i = compute_kinetic_energy(p_tr, invM); // initial kinetic energy

        // Do the calculations only for itraj=0 if isNBRA is 1
        if(isNBRA==1){
          if(itraj==0){
            hvib = ham.children[traj]->get_ham_adi();
            //hvib = projectors[traj].H() * hvib * projectors[traj];
          }
        }

        else{
          hvib = ham.children[traj]->get_ham_adi();
          //hvib = projectors[traj].H() * hvib * projectors[traj];
        }
        
        double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
        double E_f = hvib.get(new_st, new_st).real();  // final potential energy  
        double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

//        cout<<"E_i = "<<E_i<<" E_f = "<<E_f<<" T_i = "<<T_i<<" T_f = "<<T_f<<endl;

        if(T_f>=0.0){  // hop is possible - accept it
          fstates[traj] = proposed_states[traj];
        }
        else{
          fstates[traj] = initial_states[traj];
        }      
      }
      else{ 
        fstates[traj] = initial_states[traj];
      }

    }// for traj

  }// algo = 10

  else if(prms.hop_acceptance_algo==11){  // Just based on the diabatic energy levels

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];

      if(old_st != new_st){

        p_tr = p.col(traj);
        double T_i = compute_kinetic_energy(p_tr, invM); // initial kinetic energy
        double E_i = ham.children[traj]->get_ham_dia().get(old_st, old_st).real();  // initial potential energy
        double E_f = ham.children[traj]->get_ham_dia().get(new_st, new_st).real();  // final potential energy  
        double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

        if(T_f>=0.0){  // hop is possible - accept it
          fstates[traj] = proposed_states[traj];
        }
        else{
          fstates[traj] = initial_states[traj];
        }      
      }
      else{ 
        fstates[traj] = initial_states[traj];
      }

    }// for traj
  }// algo = 11


  else if(prms.hop_acceptance_algo==20){  // if rescaling momenta along the derivative coupling vector

    MATRIX dNAC(ndof, 1);

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];

      if(old_st != new_st){
        // Do the calculations only for itraj=0 if isNBRA is 1
        if(isNBRA==1){
        if(itraj==0){
        hvib = ham.children[traj]->get_ham_adi();
        //hvib = projectors[traj].H() * hvib * projectors[traj];
        }
        }
        else{
        hvib = ham.children[traj]->get_ham_adi();
        //hvib = projectors[traj].H() * hvib * projectors[traj];
        }
        double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
        double E_f = hvib.get(new_st, new_st).real();  // final potential energy  

        dNAC = 0.0;
        //for(dof = 0; dof < ndof; dof++){
        for(idof = 0; idof < ndof_active; idof++){
          dof = which_dofs[idof];
          nac = ham.children[traj]->get_dc1_adi(dof);
          //nac = projectors[traj].H() * nac * projectors[traj];

          dNAC.set(dof, 0, nac.get(old_st, new_st).real() );
        }
      
        p_tr = p.col(traj);
        if(can_rescale_along_vector(E_i, E_f, p_tr, invM, dNAC, which_dofs)){
          fstates[traj] = proposed_states[traj];
        }
        else{
          fstates[traj] = initial_states[traj];
        }
      }
      else{ 
        fstates[traj] = initial_states[traj];
      }
    }// for traj

  }// algo = 20

  else if(prms.hop_acceptance_algo==21){  // if rescaling momenta along the difference in forces

    MATRIX dF(ndof, 1);
/*
    MATRIX Fi(ndof, ntraj);
    MATRIX Ff(ndof, ntraj);
    CMATRIX _ampl_i(nst, ntraj);     // CMATRIX version of "initial_states"
    CMATRIX _ampl_f(nst, ntraj);     // CMATRIX version of "proposed_states"

    tsh_indx2vec(ham, _ampl_i, initial_states);
    Fi = ham.Ehrenfest_forces_adi(_ampl_i, 1).real();

    tsh_indx2vec(ham, _ampl_f, proposed_states);
    Ff = ham.Ehrenfest_forces_adi(_ampl_f, 1).real();
*/


    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];

      if(old_st != new_st){
        // Do the calculations only for itraj=0 if isNBRA is 1
        if(isNBRA==1){
        if(itraj==0){
        hvib = ham.children[traj]->get_ham_adi();
        //hvib = projectors[traj].H() * hvib * projectors[traj];
        }
        }
        else{
        hvib = ham.children[traj]->get_ham_adi();
        //hvib = projectors[traj].H() * hvib * projectors[traj];
        }
        double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
        double E_f = hvib.get(new_st, new_st).real();  // final potential energy        

        //for(dof = 0; dof < ndof; dof++){
        dF = 0.0;
        for(idof = 0; idof < ndof_active; idof++){
          dof = which_dofs[idof];
          df = ham.children[traj]->get_d1ham_adi(dof);
          //df = projectors[traj].H() * df * projectors[traj];
          dF.set(dof, 0, df.get(old_st, old_st).real() - df.get(new_st, new_st).real());

        }
      
        p_tr = p.col(traj);
        if(can_rescale_along_vector(E_i, E_f, p_tr, invM, dF, which_dofs)){
          fstates[traj] = proposed_states[traj];
        }
        else{
          fstates[traj] = initial_states[traj];
        }
      }
      else{ 
        fstates[traj] = initial_states[traj];
      }

    }// for traj

  }// algo = 21


  else if(prms.hop_acceptance_algo==31){  // stochastic decision based on quantum Boltzmann factors

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(itraj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy  

      double prob = boltz_factor(E_f, E_i, prms.Temperature, 1);

      double ksi = rnd.uniform(0.0, 1.0);

      if(ksi < prob ){  
        fstates[traj] = proposed_states[traj]; 
      }
      else{ 
        fstates[traj] = initial_states[traj]; 
      }

    }
  }// algo = 31

  else if(prms.hop_acceptance_algo==32){  // stochastic decision based on classical Maxwell-Boltzmann factors

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(itraj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy  

      double prob = boltz_factor(E_f, E_i, prms.Temperature, 2);

      double ksi = rnd.uniform(0.0, 1.0);

      if(ksi < prob ){  
        fstates[traj] = proposed_states[traj]; 
      }
      else{ 
        fstates[traj] = initial_states[traj]; 
      }

    }

  }// algo = 32

  else if(prms.hop_acceptance_algo==33){  // stochastic decision based on quantum probabilities

    for(itraj=0; itraj<ntraj_active; itraj++){
      traj = which_trajectories[itraj];

      int old_st = initial_states[traj];
      int new_st = proposed_states[traj];
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(itraj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy  

      double prob = boltz_factor(E_f, E_i, prms.Temperature, 3);

      double ksi = rnd.uniform(0.0, 1.0);

      if(ksi < prob ){  
        fstates[traj] = proposed_states[traj]; 
      }
      else{ 
        fstates[traj] = initial_states[traj]; 
      }

    }

  }// algo = 33


  return fstates;

}


vector<int> accept_hops(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd ){

    int ntraj = q.n_cols;
    vector<int> which_trajectories(ntraj);

    for(int i=0; i<ntraj; i++){ which_trajectories[i] = i; }

    return accept_hops(prms, q, p, invM, C, /*projectors,*/ ham, proposed_states, initial_states, rnd, which_trajectories);

}





vector<int> where_can_we_hop(int traj, dyn_control_params& prms,
       MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& act_states, Random& rnd){
/**
   This function gives a list of indices of states to which a trajectory of index `traj` can
   hop to from its current active state
   We want to determine only the non-trivial transitions (that is everything  other than hopping onto itself)

   act_states = [ state_0,  state_1, ... <state_traj>, .... ] is the array of the current state indices for each trajectory
 
   So, we are interested in the particular trajectory `traj` and we agoing to check onto which other states can that trajectory hop

   We don't care about all other values in act_states, q, p, etc. ... other than corresponding to the index `traj`

*/

    int ntraj = act_states.size();
    int nstates = Coeff.n_rows;

    vector<int> proposed_states(ntraj, -1);
    vector<int> new_states(ntraj, -1);
    vector<int> which_trajectories(1, traj);
    vector<int> all_possible_hops;

    /// Determine whereto we can hop from the current state:
    for(int istate=0; istate<nstates; istate++){

      if(istate != act_states[traj]){

        proposed_states[traj] = istate;

        /// Decide if we can accept the transitions, the function below only checks the hopping for a single trajectory `traj`
        /// other elements of the input and output vector<int> variables (old_states, new_states, proposed_states) are irrelevant
        /// the variable `which_trajectories` instructs to handle only the current trajectory
        new_states = accept_hops(prms, q, p, invM, Coeff, /*projectors,*/ ham, proposed_states, act_states, rnd, which_trajectories);

        if(new_states[traj]!=act_states[traj]){
          all_possible_hops.push_back(new_states[traj]);
        }

      }/// if not the current active state 

    }/// for istate

    return all_possible_hops;
}




void handle_hops_nuclear(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors,*/
       nHamiltonian& ham, vector<int>& new_states, vector<int>& old_states){
/**
  This function changes the nuclear dynamical variables after successful or frustrated hops

  options:
  0 - don't rescale

  100 - based on adiabatic energy, don't reverse on frustrated hops
  101 - based on adiabatic energy, reverse on frustrated hops
  110 - based on diabatic energy, don't reverse on frustrated hops
  111 - based on diabatic energy, reverse on frustrated hops

  200 - along derivative coupling vectors, don't reverse on frustrated hops
  201 - along derivative coupling vectors, reverse on frustrated hops
  210 - along difference of state-specific forces, don't reverse on frustrated hops
  211 - along difference of state-specific forces, reverse on frustrated hops

*/

  vector<int>& which_dofs = prms.quantum_dofs;

  int ndof = q.n_rows;
  int n_active_dof = which_dofs.size(); 
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, idof, dof, i;
  int isNBRA = prms.isNBRA;
  MATRIX p_tr(ndof, 1);
  CMATRIX hvib(nst, nst);
  CMATRIX nac(nst, nst);
  CMATRIX df(nst, nst);


  if(prms.momenta_rescaling_algo==0){  // Don't rescale at all
    ; ;
  }// algo = 0

  else if(prms.momenta_rescaling_algo==100 || prms.momenta_rescaling_algo==101){  // rescale momenta uniformly based on adiabatic energies

    for(traj=0; traj<ntraj; traj++){

      int old_st = old_states[traj];
      int new_st = new_states[traj];

      p_tr = p.col(traj);

      double T_i = compute_kinetic_energy(p_tr, invM, which_dofs); // initial kinetic energy
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(traj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy  
      double T_f = T_i + E_i - E_f;             // predicted final kinetic energy


      double scl_fac = 1.0;

      if(T_f>=0.0){  scl_fac = std::sqrt(T_f/T_i);   }
      else{
        if(prms.momenta_rescaling_algo==100){  scl_fac = 1.0; }
        else if(prms.momenta_rescaling_algo==101){  scl_fac = -1.0; }
      }      

      //for(dof = 0; dof < ndof; dof++){   p.scale(dof, traj, scl_fac);  }
      for(idof = 0; idof < n_active_dof; idof++){   
        dof = which_dofs[idof];
        p.scale(dof, traj, scl_fac);  
      }

    }
  }// algo = 100 || 101

  else if(prms.momenta_rescaling_algo==110 || prms.momenta_rescaling_algo==111){  // rescale momenta uniformly based on diabatic energies

    for(traj=0; traj<ntraj; traj++){

      int old_st = old_states[traj];
      int new_st = new_states[traj];

      p_tr = p.col(traj);
      double T_i = compute_kinetic_energy(p_tr, invM, which_dofs); // initial kinetic energy
      double E_i = ham.children[traj]->get_ham_dia().get(old_st, old_st).real();  // initial potential energy
      double E_f = ham.children[traj]->get_ham_dia().get(new_st, new_st).real();  // final potential energy  
      double T_f = T_i + E_i - E_f;             // predicted final kinetic energy

      double scl_fac = 1.0;

      if(T_f>=0.0){  scl_fac = std::sqrt(T_f/T_i);   }
      else{
        if(prms.momenta_rescaling_algo==110){  scl_fac = 1.0; }
        else if(prms.momenta_rescaling_algo==111){  scl_fac = -1.0; }
      }      

      //for(dof = 0; dof < ndof; dof++){   p.scale(dof, traj, scl_fac);  }
      for(idof = 0; idof < n_active_dof; idof++){   
        dof = which_dofs[idof];
        p.scale(dof, traj, scl_fac);  
      }


    }
  }// algo = 110 || algo = 111


  else if(prms.momenta_rescaling_algo==200 || prms.momenta_rescaling_algo==201){  // rescale momenta along the derivative coupling vector

    MATRIX dNAC(ndof, 1);
    int do_reverse;

    if(prms.momenta_rescaling_algo==200){ do_reverse = 0; }
    else if(prms.momenta_rescaling_algo==201){ do_reverse = 1; }

    for(traj=0; traj<ntraj; traj++){

      int old_st = old_states[traj];
      int new_st = new_states[traj];
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(traj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy  

      dNAC = 0.0;
      for(idof = 0; idof < n_active_dof; idof++){
        dof = which_dofs[idof];
        nac = ham.children[traj]->get_dc1_adi(dof);
        //nac = projectors[traj].H() * nac * projectors[traj];
        dNAC.set(dof, 0, nac.get(old_st, new_st).real() );
      }

      p_tr = p.col(traj);       
      rescale_along_vector(E_i, E_f, p_tr, invM, dNAC, do_reverse, which_dofs); 

      for(dof = 0; dof < ndof; dof++){   p.set(dof, traj, p_tr.get(dof, 0));     }


    }// for traj

  }// algo = 200 || 201

  else if(prms.momenta_rescaling_algo==210 || prms.momenta_rescaling_algo==211 ){  // rescale momenta along the difference in forces

    MATRIX dF(ndof, 1);
/*
    MATRIX Fi(ndof, ntraj);
    MATRIX Ff(ndof, ntraj);
    CMATRIX _ampl_i(nst, ntraj);     // CMATRIX version of "initial_states"
    CMATRIX _ampl_f(nst, ntraj);     // CMATRIX version of "proposed_states"


    tsh_indx2vec(ham, _ampl_i, old_states);
    Fi = ham.Ehrenfest_forces_adi(_ampl_i, 1).real();

    tsh_indx2vec(ham, _ampl_f, new_states);
    Ff = ham.Ehrenfest_forces_adi(_ampl_f, 1).real();

*/

    int do_reverse;

    if(prms.momenta_rescaling_algo==210){ do_reverse = 0; }
    else if(prms.momenta_rescaling_algo==211){ do_reverse = 1; }


    for(traj=0; traj<ntraj; traj++){

      int old_st = old_states[traj];
      int new_st = new_states[traj];
      // Do the calculations only for itraj=0 if isNBRA is 1
      if(isNBRA==1){
      if(traj==0){
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      }
      else{
      hvib = ham.children[traj]->get_ham_adi();
      //hvib = projectors[traj].H() * hvib * projectors[traj];
      }
      double E_i = hvib.get(old_st, old_st).real();  // initial potential energy
      double E_f = hvib.get(new_st, new_st).real();  // final potential energy        

      dF = 0.0;
      for(idof = 0; idof < n_active_dof; idof++){
        dof = which_dofs[idof];
        df = ham.children[traj]->get_d1ham_adi(dof);
        //df = projectors[traj].H() * df * projectors[traj];
        dF.set(dof, 0, df.get(old_st, old_st).real() - df.get(new_st, new_st).real());
      }

      p_tr = p.col(traj);       
      rescale_along_vector(E_i, E_f, p_tr, invM, dF, do_reverse, which_dofs); 

      for(dof = 0; dof < ndof; dof++){   p.set(dof, traj, p_tr.get(dof, 0));     }
      
    }// for traj

  }// algo = 210 || 211


}// handle_hops_nuclear



}// namespace libdyn
}// liblibra

