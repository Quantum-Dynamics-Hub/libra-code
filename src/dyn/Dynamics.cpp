/*********************************************************************************
* Copyright (C) 2019-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dynamics.cpp
  \brief The file implements the general framework to run:
   * adiabatic dynamics (Verlet)
   * nonadiabatic Ehrenfest dynamics
   * nonadiabatic TSH dynamics 
   * include thermostat, if needed
   * include decoherence
   * include quantum nuclear effect (ETHD) 
   * include phase corrections/state tracking in NA-MD
   * the same framework for multiple trajectories
   * include coupled-trajectories methods (planned)
   * enable the NBRA-like calculations as well as non-NBRA
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "electronic/libelectronic.h"
#include "Dynamics.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"
#include "dyn_ham.h"
#include "../calculators/NPI.h"


/// liblibra namespace
namespace liblibra{

using namespace libcalculators;

/// libdyn namespace
namespace libdyn{



namespace bp = boost::python;


void aux_get_transforms(CMATRIX** Uprev, nHamiltonian& ham){

  // For adiabatic representation only:
  // Save the previous orbitals info - in case we need to
  // do either phase correction of state tracking

  int ntraj = ham.children.size();

  for(int traj=0; traj<ntraj; traj++){
    *Uprev[traj] = ham.children[traj]->get_basis_transform();  
  }

}




//vector<CMATRIX> compute_St(nHamiltonian& ham, CMATRIX** Uprev){
//vector<CMATRIX> compute_St(nHamiltonian& ham, vector<CMATRIX>& Uprev, int isNBRA){
vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev, int isNBRA){

  return compute_St(&ham, &ham_prev, isNBRA);
}

vector<CMATRIX> compute_St(nHamiltonian* ham, nHamiltonian* ham_prev, int isNBRA){
/**
  This function computes the time-overlap matrices for all trajectories

*/

  int nst = ham->nadi;
  int ntraj = ham->children.size();

  vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));

  if(isNBRA==1){
    St[0] = ham_prev->children[0]->get_basis_transform().H() * ham->children[0]->get_basis_transform();    
    // AVA: temporarily comment on 11/25/2022
    ham->children[0]->set_time_overlap_adi_by_val(St[0]);
  }
  else{
    for(int traj=0; traj<ntraj; traj++){
      St[traj] = ham_prev->children[traj]->get_basis_transform().H() * ham->children[traj]->get_basis_transform();
      ham->children[traj]->set_time_overlap_adi_by_val(St[traj]);      
    }
  }
  return St;

}

//vector<CMATRIX> compute_St(nHamiltonian& ham, vector<CMATRIX>& Uprev){
vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev){

  int is_nbra = 0;
  return compute_St(&ham, &ham_prev, is_nbra);

}

vector<CMATRIX> compute_St(nHamiltonian* ham, int isNBRA){
/**
  This function computes the time-overlap matrices for all trajectories
*/

  int nst = ham->nadi;
  int ntraj = ham->children.size();

  vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));
  if(isNBRA==1){
    St[0] = ham->children[0]->get_time_overlap_adi();
  }
  else{
    for(int traj=0; traj<ntraj; traj++){
      St[traj] = ham->children[traj]->get_time_overlap_adi();
    }
  }
  return St;

}


vector<CMATRIX> compute_St(nHamiltonian& ham, int isNBRA){

  return compute_St(&ham, isNBRA);
}

vector<CMATRIX> compute_St(nHamiltonian& ham){
  int is_nbra = 0;

  return compute_St(&ham, is_nbra);
}



void apply_afssh(dyn_variables& dyn_var, CMATRIX& C, vector<int>& act_states, MATRIX& invM,
                nHamiltonian& ham, bp::dict& dyn_params, Random& rnd  ){

  //cout<<"In apply_afssh\n";

  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int i,j;
  int ndof = invM.n_rows;
  int nst = C.n_rows;    
  int ntraj = C.n_cols;
  int traj, dof, idof;
  int num_el = prms.num_electronic_substeps;
  double dt_el = prms.dt / num_el;

  // A-FSSH

    CMATRIX hvib_curr(nst, nst);
    CMATRIX force_full(nst, nst);
    CMATRIX force_diag(nst, nst);
    CMATRIX c_traj(nst, 1);
    CMATRIX dR_afssh(nst, nst);
    CMATRIX dP_afssh(nst, nst);



    //cout<<"Propagating moments...\n";
    //=========================== Propagate moments ===============
    for(traj=0; traj<ntraj; traj++){

      hvib_curr = ham.children[traj]->get_hvib_adi();
      c_traj = C.col(traj);

      double gamma_reset = 0.0;

      for(idof=0; idof<ndof; idof++){  

         force_full = -1.0 * ham.children[traj]->get_d1ham_adi(idof);

         for(i=0;i<nst;i++){  force_diag.set(i, i,  force_full.get(i,i) ); } 

         dR_afssh = *dyn_var.dR[traj][idof];
         dP_afssh = *dyn_var.dP[traj][idof];

         integrate_afssh_moments(dR_afssh, dP_afssh, hvib_curr, force_diag, 
                                 c_traj, 1.0/invM.get(idof,0), act_states[traj], dt_el, num_el);

         *dyn_var.dR[traj][idof] = dR_afssh; 
         *dyn_var.dP[traj][idof] = dP_afssh;         

      }// for idof
    }// for traj


    //cout<<"Computing reset and collapse probabilities...\n";

    //======================== Compute reset and collapse probabilities =========

    MATRIX gamma_reset(nst, ntraj);
    MATRIX gamma_collapse(nst, ntraj);


    for(traj=0; traj<ntraj; traj++){
      for(i=0;i<nst;i++){

        double gamma_reset_i = 0.0;
        double gamma_collapse_i = 0.0;

        for(idof=0; idof<ndof; idof++){  
        
          double dx_ii = dR_afssh.get(i, i).real();
          int as = act_states[traj];
          double f_i   = -ham.children[traj]->get_d1ham_adi(idof).get(i, i).real();
          double f_as = -ham.children[traj]->get_d1ham_adi(idof).get(as,as).real();

          gamma_reset_i -= 0.5*(f_i - f_as) * dx_ii;

          double f_ji   = -ham.children[traj]->get_d1ham_adi(idof).get(as, i).real();
          gamma_collapse_i += f_ji * dx_ii;

        }// for idof

        gamma_reset.set(i, traj, gamma_reset_i * prms.dt);
        gamma_collapse.set(i, traj, (gamma_reset_i -  2.0*fabs(gamma_collapse_i)) * prms.dt );

      }// for nst
    }// for traj


    //cout<<"Doing the collapses and resets...\n";
    //======================== Do the collapse and resets =======================

    complex<double> zero(0.0, 0.0);

    for(traj=0; traj<ntraj; traj++){

      for(i=0;i<nst;i++){

        if(i!=act_states[traj]){


          // Reset
          double ksi = rnd.uniform(0.0, 1.0);

          if(ksi < gamma_reset.get(i, traj)){
            for(idof=0;idof<ndof;idof++){
              dyn_var.dR[traj][idof]->scale(-1, i, zero);
              dyn_var.dR[traj][idof]->scale(i, -1, zero);    
              dyn_var.dP[traj][idof]->scale(-1, i, zero);
              dyn_var.dP[traj][idof]->scale(i, -1, zero);
            }// for j
          }// if ksi < gamma_reset


          // Collapse
          ksi = rnd.uniform(0.0, 1.0);
          if(ksi < gamma_collapse.get(i, traj)){
            collapse(C, traj, act_states[traj], prms.collapse_option);
          }// if ksi < gamma_collapse
          

        }// all non-active states

      }// for i
    }// for traj

    // cout<<"Done\n";
}




void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
              vector<int>& act_states,              
              nHamiltonian& ham, bp::object py_funct, bp::dict params, bp::dict dyn_params, Random& rnd){
/**
  This is a version to maintain the backward-compatibility
*/ 

  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int ntraj = q.n_cols;
  vector<Thermostat> therm(ntraj, Thermostat(prms.thermostat_params));

  compute_dynamics(q, p, invM, C, projectors, act_states, ham, py_funct, params, dyn_params, rnd, therm);

}

void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
              vector<int>& act_states,              
              nHamiltonian& ham, bp::object py_funct, bp::dict& params, bp::dict& dyn_params, Random& rnd,
              vector<Thermostat>& therm){

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    

  dyn_variables dyn_var(nst, nst, ndof, ntraj);
  compute_dynamics(q, p, invM, C, projectors, act_states, ham, py_funct, params, dyn_params, rnd, therm, dyn_var);

}


void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
              vector<int>& act_states,              
              nHamiltonian& ham, bp::object py_funct, bp::dict& params, bp::dict& dyn_params, Random& rnd,
              vector<Thermostat>& therm, dyn_variables& dyn_var){

/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x ntraj]  or [ndia x ntraj] matrix containing the electronic coordinates. The amplitudes
   are assumed to be dynamically-consistent
  \param[in,out] projectors [ntraj CMATRIX(nadi, nadi)] - the projector matrices that account for the state tracking and 
  phase correction. These matrices should be considered as the dynamical varibles, similar to quantum amplitudes. Except
  their evolution does not necessarily follow from some equations of motion, but rather from various ad hoc schemes.
  \param[in,out] act_states - vector of ntraj indices of the physical states in which each of the trajectories
  initially is (active states). 
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) 
  - its internal variables (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] params1 The Python dictionary containing the control parameters passed to this function
  \param[in] rnd The Random number generator object

  Return: propagates C, q, p and updates state variables

*/


/*

  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int i,j;
  int cdof;
  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int ntraj1;
  int nst = C.n_rows;    
  int traj, dof, idof;
  int n_therm_dofs;
  int num_el = prms.num_electronic_substeps;
  double dt_el = prms.dt / num_el;
  vector< vector<int> > perms;

  if(prms.rep_tdse==0){  dyn_var.ampl_dia = &C;  }
  else if(prms.rep_tdse==1){  dyn_var.ampl_adi = &C;  }
  dyn_var.act_states = act_states;


  MATRIX coherence_time(nst, ntraj); // for DISH
  MATRIX coherence_interval(nst, ntraj); // for DISH
  vector<int> project_out_states(ntraj); // for DISH

  vector<CMATRIX> insta_proj(ntraj, CMATRIX(nst, nst));
 
  vector<CMATRIX> Uprev;
  // Defining ntraj1 as a reference for making these matrices
  // ntraj is defined as q.n_cols as since it would be large in NBRA
  // we can define another variable like ntraj1 and build the matrices based on that.
  // We can make some changes where q is generated but this seems to be a bit easier
  if(prms.isNBRA==1){   ntraj1 = 1;  }
  else{  ntraj1 = ntraj;  }

  // Defining matrices based on ntraj1
  vector<CMATRIX> St(ntraj1, CMATRIX(nst, nst));
  vector<CMATRIX> Eadi(ntraj1, CMATRIX(nst, nst));  
  vector<MATRIX> decoherence_rates(ntraj1, MATRIX(nst, nst)); 
  vector<double> Ekin(ntraj1, 0.0);  
  vector<MATRIX> prev_ham_dia(ntraj1, MATRIX(nst, nst));
  MATRIX gamma(ndof, ntraj);
  MATRIX p_traj(ndof, 1);
  vector<int> t1(ndof, 0); for(dof=0;dof<ndof;dof++){  t1[dof] = dof; }
  vector<int> t2(1,0);
  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }
  CMATRIX c_tmp(nst, 1);

  //if(prms.decoherence_algo==2){   dyn_var.allocate_afssh(); }

  //============ Sanity checks ==================
  if(prms.ensemble==1){  
    n_therm_dofs = therm[0].Nf_t + therm[0].Nf_r;
    if(n_therm_dofs != prms.thermostat_dofs.size()){
      cout<<"Error in compute_dynamics: The number of thermostat DOFs ( currently "<<n_therm_dofs<<") must be \
      equal to the number of thermostat dofs set up by the `thermostat_dofs` parameter ( currently "
      <<prms.thermostat_dofs.size()<<")\nExiting...\n";
      exit(0);
    }
  }


 if(prms.tsh_method == 3){ // DISH
  //    prev_ham_dia[0] = ham.children[0]->get_ham_dia().real();
  //}
  //else{
    for(traj=0; traj<ntraj1; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }
  }
  //============ Update the Hamiltonian object =============
  // In case, we may need phase correction & state reordering
  // prepare the temporary files
  if(prms.rep_tdse==1){      
    if(prms.do_phase_correction || prms.state_tracking_algo > 0){

      // On-the-fly calculations, from the wavefunctions
      if(prms.time_overlap_method==0){
        Uprev = vector<CMATRIX>(ntraj, CMATRIX(nst, nst));

        for(traj=0; traj<ntraj; traj++){
          Uprev[traj] = ham.children[traj]->get_basis_transform();  
        }
      }      
    }// do_phase_correction || state_tracking_algo > 0
  }// rep == 1

  // In case, we select to compute scalar NACs from time-overlaps
  update_nacs(prms, ham);


  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  // Adding the prms.isNBRA to the propagate electronic
  update_Hamiltonian_p(prms, ham, p, invM);
  for(i=0; i<num_el; i++){
    propagate_electronic(0.5*dt_el, C, projectors, ham.children, prms.rep_tdse, prms.isNBRA);   
  }


  //============== Nuclear propagation ===================
  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        p.scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }

  p = p + aux_get_forces(prms, dyn_var, ham) * 0.5 * prms.dt;

  // Kinetic constraint
  for(cdof = 0; cdof < prms.constrained_dofs.size(); cdof++){   
    p.scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }



  if(prms.entanglement_opt==22){
    gamma = ETHD3_friction(q, p, invM, prms.ETHD3_alpha, prms.ETHD3_beta);
  }
  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * prms.dt ); 

      if(prms.entanglement_opt==22){
        q.add(dof, traj,  invM.get(dof,0) * gamma.get(dof,traj) * prms.dt ); 
      }

    }
  }


  // Recompute the matrices at the new geometry and apply any necessary fixes 
  update_Hamiltonian_q(prms, q, ham, py_funct, params);
  update_Hamiltonian_q_ethd(prms, q, p, ham, py_funct, params, invM);


  // Apply phase correction and state reordering as needed
  if(prms.rep_tdse==1){

    if(prms.state_tracking_algo > 0 || prms.do_phase_correction){

      /// Compute the time-overlap directly, using previous MO vectors
      if(prms.time_overlap_method==0){    
        St = compute_St(ham, Uprev, prms.isNBRA);
      }
      /// Read the existing time-overlap
      else if(prms.time_overlap_method==1){    St = compute_St(ham, prms.isNBRA);      }
      Eadi = get_Eadi(ham);           // these are raw properties
      //update_projectors(prms, projectors, Eadi, St, rnd);

      perms = compute_permutations(prms, Eadi, St, rnd);
      insta_proj = compute_projectors(prms, St, perms);
//      insta_proj = compute_projectors(prms, Eadi, St, rnd);
      
      /// Adiabatic Amplitudes
      for(traj=0; traj<ntraj; traj++){
        t2[0] = traj;
        pop_submatrix(C, c_tmp, t3, t2);
        c_tmp = insta_proj[traj] * c_tmp;
        push_submatrix(C, c_tmp, t3, t2);
      }

      /// Adiabatic states are permuted
      act_states = permute_states(perms, act_states);

    }
  }// rep_tdse == 1


  // In case, we select to compute scalar NACs from time-overlaps
  update_nacs(prms, ham);

  // NVT dynamics
  if(prms.ensemble==1){  
    for(traj=0; traj<ntraj; traj++){
      t2[0] = traj; 
      pop_submatrix(p, p_traj, t1, t2);
      double ekin = compute_kinetic_energy(p_traj, invM, prms.thermostat_dofs);
      therm[traj].propagate_nhc(prms.dt, ekin, 0.0, 0.0);
    }

  }

  p = p + aux_get_forces(prms, dyn_var, ham) * 0.5 * prms.dt;


  // Kinetic constraint
  for(cdof=0; cdof<prms.constrained_dofs.size(); cdof++){   
    p.scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }

  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        p.scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }

  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  update_nacs(prms, ham);
  update_Hamiltonian_p(prms, ham, p, invM);
  for(i=0; i<num_el; i++){
    propagate_electronic(0.5*dt_el, C, projectors, ham.children, prms.rep_tdse, prms.isNBRA);
  }

  CMATRIX Hvib(ham.children[0]->nadi, ham.children[0]->nadi);  
  Hvib = ham.children[0]->get_hvib_adi();


  //============== Begin the TSH part ===================

  // To be able to compute transition probabilities, compute the corresponding amplitudes
  // This transformation is between diabatic and raw adiabatic representations
  CMATRIX Coeff(nst, ntraj);   // this is adiabatic amplitudes
  
  Coeff = C;


  //================= Update decoherence rates & times ================
  //MATRIX decoh_rates(*prms.decoh_rates);
//    exit(0);
  if(prms.decoherence_times_type==-1){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = 0.0;   }
  }

  /// mSDM
  /// Just use the plain times given from the input, usually the
  /// mSDM formalism
  else if(prms.decoherence_times_type==0){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = *prms.decoherence_rates;   }
  }

  /// Compute the dephasing rates according the original energy-based formalism
  else if(prms.decoherence_times_type==1){
    Eadi = get_Eadi(ham); 
    Ekin = compute_kinetic_energies(p, invM);
    decoherence_rates = edc_rates(Eadi, Ekin, prms.decoherence_C_param, prms.decoherence_eps_param, prms.isNBRA);       
  }

  else if(prms.decoherence_times_type==2){
    decoherence_rates = schwartz_1(prms, Coeff, projectors, ham, *prms.schwartz_decoherence_inv_alpha); 
  }

  else if(prms.decoherence_times_type==3){
    decoherence_rates = schwartz_2(prms, projectors, ham, *prms.schwartz_decoherence_inv_alpha); 
  }


  ///== Optionally, apply the dephasing-informed correction ==
  if(prms.dephasing_informed==1){
    Eadi = get_Eadi(ham); 
    MATRIX ave_gaps(*prms.ave_gaps);
    dephasing_informed_correction(decoherence_rates, Eadi, ave_gaps, prms.isNBRA);
  }

  //============ Apply decoherence corrections ==================
  // SDM and alike methods
  if(prms.decoherence_algo==0){
    Coeff = sdm(Coeff, prms.dt, act_states, decoherence_rates, prms.sdm_norm_tolerance, prms.isNBRA);
  }
  // BCSH
  else if(prms.decoherence_algo==3){ 
    *dyn_var.reversal_events = wp_reversal_events(p, invM, act_states, ham, projectors, prms.dt);
    Coeff = bcsh(Coeff, prms.dt, act_states, *dyn_var.reversal_events);
  }

  else if(prms.decoherence_algo==4){
//    MATRIX decoh_rates(ndof, ntraj);
    Coeff = mfsd(p, Coeff, invM, prms.dt, decoherence_rates, ham, rnd, prms.isNBRA);
  }

//  exit(0);

  //========= Use the resulting amplitudes to do the hopping =======
  // Adiabatic dynamics
  if(prms.tsh_method==-1){ ;; } 

  // FSSH, GFSH, MSSH
  else if(prms.tsh_method == 0 || prms.tsh_method == 1 || prms.tsh_method == 2){

    // Compute hopping probabilities
    vector<MATRIX> g( hop_proposal_probabilities(prms, q, p, invM, Coeff, projectors, ham, prev_ham_dia) );

    // Propose new discrete states    
    vector<int> prop_states( propose_hops(g, act_states, rnd) );


    // Decide if to accept the transitions (and then which)
    vector<int> old_states(act_states);
    act_states = accept_hops(prms, q, p, invM, Coeff, projectors, ham, prop_states, act_states, rnd);    

    // Velocity rescaling
    handle_hops_nuclear(prms, q, p, invM, Coeff, projectors, ham, act_states, old_states);



    if(prms.decoherence_algo==1){
      // Instantaneous decoherence 
      instantaneous_decoherence(Coeff, act_states, prop_states, old_states, 
          prms.instantaneous_decoherence_variant, prms.collapse_option);
    } 
    else if(prms.decoherence_algo==2){
//      exit(0);
      apply_afssh(dyn_var, C, act_states, invM, ham, dyn_params, rnd);

    }// AFSSH
        
  }// tsh_method == 0, 1, 2, 4
  
  // DISH
  else if(prms.tsh_method == 3 ){

    /// Advance coherence times
    coherence_time.add(-1, -1, prms.dt);

    /// Older version:

    /// New version, as of 8/3/2020
    vector<int> old_states(act_states);
    act_states = dish(prms, q, p, invM, Coeff, projectors, ham, act_states, coherence_time, decoherence_rates, rnd);

    /// Velocity rescaling
    handle_hops_nuclear(prms, q, p, invM, Coeff, projectors, ham, act_states, old_states);


  }// tsh_method == 3

  else{   cout<<"tsh_method == "<<prms.tsh_method<<" is undefined.\nExiting...\n"; exit(0);  }


  // Need the inverse
  //Coeff = transform_amplitudes(prms.rep_tdse, prms.rep_sh, C, ham);
  C = Coeff;


  project_out_states.clear();

  if(prms.rep_tdse==1){      
    if(prms.do_phase_correction || prms.state_tracking_algo > 0){
      if(prms.time_overlap_method==0){  Uprev.clear();  }
    }
  }

//    exit(0);
  St.clear();
  Eadi.clear();
  decoherence_rates.clear();
  Ekin.clear();
  prev_ham_dia.clear();

  dyn_var.act_states = act_states;

*/

}

void apply_projectors(CMATRIX& C, vector<CMATRIX>& proj){

  int ntraj = proj.size();
  int nst = C.n_rows;

  int i;
  vector<int> t2(1,0);
  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }
  CMATRIX c_tmp(nst, 1);

  for(int traj=0; traj<ntraj; traj++){
    t2[0] = traj;
    pop_submatrix(C, c_tmp, t3, t2);
    c_tmp = proj[traj] * c_tmp;
    push_submatrix(C, c_tmp, t3, t2);
  }
}

void propagate_electronic(double dt, CMATRIX& C, vector<CMATRIX*>& proj, nHamiltonian& ham, nHamiltonian& ham_prev, dyn_control_params& prms, Random& rnd){

  propagate_electronic(dt, C, proj, &ham, &ham_prev, prms, rnd);

}


void propagate_electronic(double dt, CMATRIX& Coeff, vector<CMATRIX*>& proj, nHamiltonian* Ham, nHamiltonian* Ham_prev, dyn_control_params& prms, Random& rnd){

  int itraj;
  int ntraj = Ham->children.size();
  int nst = Coeff.n_rows;
  CMATRIX C(nst, 1);

  int rep = prms.rep_tdse;
  int method = prms.electronic_integrator;


//  cout<<"ntraj = "<<ntraj<<" nst = "<<nst<<" rep = "<<rep<<" method = "<<method<<endl;
//  exit(0);
//  int i;
//  vector<int> t2(1,0);
//  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }
//  CMATRIX c_tmp(nst, 1);

  //==================== Determine the state-consistency transformations =================
  // Apply phase correction and state reordering as needed
/*
  vector<CMATRIX> insta_proj(ntraj, CMATRIX(nst, nst));
  vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));
  vector<CMATRIX> T_new(ntraj, CMATRIX(nst, nst));
*/
//  cout<<"Point 0\n"; exit(0);
//  if(prms.state_tracking_algo > 0 || prms.do_phase_correction){

    ///==================== Time-overlaps =====================
    //vector<CMATRIX> St(ntraj, CMATRIX(nst, nst));
    /// Compute the time-overlap directly, using previous MO vectors
//    if(prms.time_overlap_method==0){  St = compute_St(Ham, Ham_prev, prms.isNBRA);   }
    /// Read the existing time-overlap
//    else if(prms.time_overlap_method==1){    St = compute_St(Ham, prms.isNBRA);   }

//    cout<<"Point 1\n"; exit(0);

    ///==================== Energies ===========================
    //vector<CMATRIX> Eadi(ntraj, CMATRIX(nst, nst));
    //Eadi = get_Eadi(Ham);           // these are raw properties

    ///==== Compute permutations and insta projections =========
    //vector< vector<int> > perms;
    //perms = compute_permutations(prms, Eadi, St, rnd);

    //vector<CMATRIX> insta_proj(ntraj, CMATRIX(nst, nst));
    //insta_proj = compute_projectors(prms, St, perms);

    ///New approach to compute projectors:
    //CMATRIX TP(nst, nst);
   
    //for(itraj=0; itraj<ntraj;itraj++){        
    //    TP = (*proj[itraj]).H() * St[itraj];
    //    FullPivLU_inverse(TP, T_new[itraj]);
    //}
    

    //if(prms.rep_tdse==1){   apply_projectors(C, insta_proj); }
    /// Adiabatic states are permuted
//    act_states = permute_states(perms, act_states);
//  }

//  cout<<"Here\n"; exit(0);

  ///======================== Now do the integration of the TD-SE ===================
  for(itraj=0; itraj<ntraj; itraj++){

    C = Coeff.col(itraj);

    int traj1 = itraj;  if(method >=100 and method <200){ traj1 = 0; }

    nHamiltonian* ham = Ham->children[traj1];
    nHamiltonian* ham_prev = Ham_prev->children[traj1];


    CMATRIX P(ham->nadi, ham->nadi);
    proj[itraj]->load_identity();
    CMATRIX& T = *proj[itraj];
    CMATRIX T_new(ham->nadi, ham->nadi);

    P = ham->get_time_overlap_adi();  // U_old.H() * U;
    P = orthogonalized_T( P );
    FullPivLU_inverse(P, T_new);
    


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

    if(method==-1 || method==100){
/*
      CMATRIX Hvib(ham->nadi, ham->nadi);  
      CMATRIX T(nst, nst); T = *proj[itraj];
      CMATRIX iT(nst, nst);
      CMATRIX TP(nst, nst);
      CMATRIX T_new(nst, nst); 
      CMATRIX iTP(nst, nst);

      TP = T.H() * St[itraj];
      FullPivLU_inverse(TP, T_new);

      //Hvib = T.H() * ham_prev->get_hvib_adi() * T;
      //Hvib += T_new.H() * ham->get_hvib_adi() * T_new;
      Hvib = T.H() * ham_prev->get_ham_adi() * T;
      Hvib += T_new.H() * ham->get_ham_adi() * T_new;
      Hvib *= 0.5;
  
      FullPivLU_inverse(T, iT);
      C = iT * C;

      //cout<<"Integration with projector T = "; T.show_matrix();
      propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs
      //propagate_electronic_nonHermitian(dt, C, Hvib);
      //CMATRIX I(nst, nst); I.load_identity();
      //propagate_electronic(dt, C, Hvib, I);

      C = T_new * C;

      *proj[itraj] = T_new;
*/
    }
 
    else if(method==0 || method==100){

//      cout<<"============= Start integration =============\n";
      // A new approach
      CMATRIX Hvib(ham->nadi, ham->nadi);
      CMATRIX A(ham->nadi, ham->nadi);
      CMATRIX B(ham->nadi, ham->nadi);



//      TP = T.H() * P;
//      FullPivLU_inverse(TP, T_new);
//       T_new = iP;

//      cout<<"T_new = \n"; T_new.show_matrix();

/*
      CMATRIX st2(ham->nadi, ham->nadi); st2 = T_new.H() * T_new;
      CMATRIX st2_half(ham->nadi, ham->nadi);
      CMATRIX st2_i_half(ham->nadi, ham->nadi);
      sqrt_matrix(st2, st2_half, st2_i_half);
      T_new = T_new * st2_i_half;
*/
//      T_new = orthogonalized_T( T_new );

//      *proj[itraj] = T_new; // *proj[itraj] * P; // load_identity();
//      cout<<"T(t+dt) = \n"; T_new.show_matrix();

//      cout<<"In electronic evolution \n"; T_new.show_matrix();

      //Hvib = T.H() * ham_prev->get_hvib_adi() * T;
      //Hvib += T_new.H() * ham->get_hvib_adi() * T_new;
      //Hvib = T.H() * ham_prev->get_ham_adi() * T;
      //Hvib += T_new.H() * ham->get_ham_adi() * T_new;


//      CMATRIX u(ham->nadi, ham->nadi);
//      CMATRIX I(ham->nadi, ham->nadi); I.load_identity();
      
      Hvib = T.H() * ham_prev->get_ham_adi() * T;
      //Hvib = ham_prev->get_ham_adi();
      A = exp_(Hvib, complex<double>(0.0, -0.5*dt) );

//      cout<<"H_prev = \n"; Hvib.show_matrix();

      Hvib = T_new.H() * ham->get_ham_adi() * T_new;
      //Hvib = ham->get_ham_adi();
//      cout<<"H = \n"; Hvib.show_matrix();

      B = exp_(Hvib, complex<double>(0.0, -0.5*dt) );

/*
      u = I - 0.5*u;
      CMATRIX iu(ham->nadi, ham->nadi);
      FullPivLU_inverse(u, iu);
*/


//      C = P.H() * A * iP.H() * B * iP * A * C;
//      C = B * iP * A * C;
//      cout<<"C(t) = \n"; C.show_matrix();
//        C = B * iP * A * C;
       C = B * T_new * A * C;

//        C = B * A * C;
//      C = A * C;

//      cout<<"C(t+dt) = \n"; C.show_matrix();
//      C = T_new * 

      // 
//      *proj[itraj] = *proj[itraj] * P; // load_identity();

//      cout<<"============= End of integration =============\n";

    }
    else if(method==1 || method==101){
      CMATRIX Hvib(ham->nadi, ham->nadi);
      Hvib = ham_prev->get_hvib_adi();
      //Hvib += insta_proj[itraj] * ham->get_hvib_adi() * insta_proj[itraj].H();
      Hvib *= 0.5;

//      cout<<"Integration with Hvib = "; Hvib.show_matrix();
      propagate_electronic_rot(dt, C, Hvib);  // in this case C - adiabatic coeffs

      //apply_projectors(C, insta_proj);
      //C = insta_proj[itraj] * C;

      //Hvib = ham->get_hvib_adi();

      //propagate_electronic_rot(0.5*dt, C, Hvib);

    }



    else if(method==2 || method==102){
      // Local diabatization
      CMATRIX U_old(ham->nadi, ham->nadi);   U_old = ham_prev->get_basis_transform();
      CMATRIX U(ham->nadi, ham->nadi);   U = ham->get_basis_transform();

      CMATRIX st(ham->nadi, ham->nadi); st = U_old.H() * U;
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

    }// method == 2

  }// rep == 1 - adiabatic


  *proj[itraj] = T_new;

  // Insert the propagated result back
  for(int st=0; st<nst; st++){  Coeff.set(st, itraj, C.get(st, 0));  }

  }// for itraj - all trajectories


}

/*
void propagate_electronic(double dt, CMATRIX& C, vector<CMATRIX*>& proj, vector<nHamiltonian*>& ham, vector<nHamiltonian*>& ham_prev, dyn_control_params& prms, Random& rnd){

  int rep = prms.rep_tdse;
  int method = prms.electronic_integrator;


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
    propagate_electronic(dt, ctmp, proj, ham[traj1], ham_prev[traj1], prms, rnd);

    // Insert the propagated result back
    for(int st=0; st<nst; st++){  C.set(st, traj, ctmp.get(st, 0));  }
  }


}
*/

/*
vector<int> update_active_states(vector<int>& act_states, vector<CMATRIX*>& T){

//   act_states[itraj] - index of the active adiabatic state for the trajectory `itraj`
//                       in terms of the time-ordered states


  int ntraj = act_state.size();
  int nst = T[0]->n_cols;
  vector<int> res( act_states);
  
  CMATRIX t2(nst,nst);
  CMATRIX t2_half(nst, nst);
  CMATRIX t2_i_half(nst, nst);
  CMATRIX coeff(nst, 1);
  complex<double> max_val;  
 
  for(int itraj=0; itraj<ntraj; itraj++){
    CMATRIX& t = *T[itraj]; 
    t2 = t.H() * t; 
    sqrt_matrix(t2, t2_half, t2_i_half);
    t2 = t * t2_i_half; 
    coeff.set(act_state[itraj], 0, complex<double>(1.0, 0.0));
    coeff = t2 * coeff;
    coeff.max_col_elt(0, max_val, res[itraj]);
  }// for itraj

  return res;
}
*/

void compute_dynamics(dyn_variables& dyn_var, bp::dict dyn_params,
              nHamiltonian& ham, nHamiltonian& ham_aux, bp::object py_funct, bp::dict params,  Random& rnd,
              vector<Thermostat>& therm){
/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x ntraj]  or [ndia x ntraj] matrix containing the electronic coordinates. The amplitudes
   are assumed to be dynamically-consistent
  \param[in,out] projectors [ntraj CMATRIX(nadi, nadi)] - the projector matrices that account for the state tracking and 
  phase correction. These matrices should be considered as the dynamical varibles, similar to quantum amplitudes. Except
  their evolution does not necessarily follow from some equations of motion, but rather from various ad hoc schemes.
  \param[in,out] act_states - vector of ntraj indices of the physical states in which each of the trajectories
  initially is (active states). 
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) 
  - its internal variables (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] params1 The Python dictionary containing the control parameters passed to this function
  \param[in] rnd The Random number generator object

  Return: propagates C, q, p and updates state variables

*/

  //========= Control parameters variables ===========
  dyn_control_params prms;
  prms.set_parameters(dyn_params);

  int num_el = prms.num_electronic_substeps;
  double dt_el = prms.dt / num_el;

  //======= Parameters of the dyn variables ==========
  int ndof = dyn_var.ndof; 
  int ntraj = dyn_var.ntraj;   
  int nadi = dyn_var.nadi;
  int ndia = dyn_var.ndia;

  //cout<<"Start dyn func\n"; exit(0);

  int nst;
  if(prms.rep_tdse==0){ nst = ndia; }
  else if(prms.rep_tdse==1){ nst = nadi; }

  //if(prms.decoherence_algo==2){   dyn_var.allocate_afssh(); }

  //========== Aliases ===============================
  CMATRIX Cact(nst, nst); 
  
//  CMATRIX& Cadi = dyn_var.ampl_adi;
//  CMATRIX& Cdia = *dyn_var.ampl_dia;
//  CMATRIX Cact; // = *dyn_var.ampl_dia;  // active rep amplitudes
  if(prms.rep_tdse==0){ Cact  = *dyn_var.ampl_dia; } //*dyn_var.ampl_dia; }
  else if(prms.rep_tdse==1){ Cact = *dyn_var.ampl_adi; } //*dyn_var.ampl_adi; }

//  cout<<"Cadi = \n"; Cadi.show_matrix(); 
//  cout<<"dyn_var_adi = \n"; dyn_var.ampl_adi->show_matrix();
//  cout<<"Cdia = \n"; Cdia.show_matrix();
//  cout<<"dyn_var_dia = \n"; dyn_var.ampl_dia->show_matrix();
//  cout<<"Cact = \n"; Cact.show_matrix();

  vector<int>& act_states = dyn_var.act_states;
  MATRIX& q = *dyn_var.q;
  MATRIX& p = *dyn_var.p;
  MATRIX& invM = *dyn_var.iM;
  
  //nHamiltonian ham_prev(ham);

  //======== General variables =======================
  int i,j, cdof, traj, dof, idof, ntraj1, n_therm_dofs;
//  vector< vector<int> > perms;


  MATRIX coherence_time(nst, ntraj);     // for DISH
  MATRIX coherence_interval(nst, ntraj); // for DISH
  vector<int> project_out_states(ntraj);  // for DISH

//  vector<CMATRIX> insta_proj(ntraj, CMATRIX(nst, nst));
 
  //vector<CMATRIX> Uprev;
  // Defining ntraj1 as a reference for making these matrices
  // ntraj is defined as q.n_cols as since it would be large in NBRA
  // we can define another variable like ntraj1 and build the matrices based on that.
  // We can make some changes where q is generated but this seems to be a bit easier
  if(prms.isNBRA==1){   ntraj1 = 1;  }
  else{  ntraj1 = ntraj;  }

  // Defining matrices based on ntraj1
//  vector<CMATRIX> St(ntraj1, CMATRIX(nst, nst));
  vector<CMATRIX> Eadi(ntraj1, CMATRIX(nst, nst));  
  vector<MATRIX> decoherence_rates(ntraj1, MATRIX(nst, nst)); 
  vector<double> Ekin(ntraj1, 0.0);  
  vector<MATRIX> prev_ham_dia(ntraj1, MATRIX(nst, nst));
  MATRIX gamma(ndof, ntraj);
  MATRIX p_traj(ndof, 1);
  vector<int> t1(ndof, 0); for(dof=0;dof<ndof;dof++){  t1[dof] = dof; }
  vector<int> t2(1,0);
  vector<int> t3(nst, 0); for(i=0;i<nst;i++){  t3[i] = i; }
//  CMATRIX c_tmp(nst, 1);
  //MATRIX F_eff(ndof, ntraj);


  //============ Sanity checks ==================
  if(prms.ensemble==1){  
    n_therm_dofs = therm[0].Nf_t + therm[0].Nf_r;
    if(n_therm_dofs != prms.thermostat_dofs.size()){
      cout<<"Error in compute_dynamics: The number of thermostat DOFs ( currently "<<n_therm_dofs<<") must be \
      equal to the number of thermostat dofs set up by the `thermostat_dofs` parameter ( currently "
      <<prms.thermostat_dofs.size()<<")\nExiting...\n";
      exit(0);
    }
  }



 if(prms.tsh_method == 3){ // DISH
  //    prev_ham_dia[0] = ham.children[0]->get_ham_dia().real();
  //}
  //else{
    for(traj=0; traj<ntraj1; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }
  }




  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  // Adding the prms.isNBRA to the propagate electronic

  //============== Nuclear propagation ===================
  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        p.scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }

  //F_eff = aux_get_forces(prms, dyn_var, ham); 
  //p = p + F_eff * 0.5 * prms.dt;
  *dyn_var.p = *dyn_var.p + 0.5*prms.dt* (*dyn_var.f);

  // Kinetic constraint
  for(cdof = 0; cdof < prms.constrained_dofs.size(); cdof++){   
    p.scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }



  if(prms.entanglement_opt==22){
    gamma = ETHD3_friction(q, p, invM, prms.ETHD3_alpha, prms.ETHD3_beta);
  }
  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * prms.dt ); 

      if(prms.entanglement_opt==22){
        q.add(dof, traj,  invM.get(dof,0) * gamma.get(dof,traj) * prms.dt ); 
      }

    }
  }


  // Recompute the matrices at the new geometry and apply any necessary fixes 
  ham_aux.copy_content(ham);
  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 0);
  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1);

//  cout<<"Exit here\n"; exit(0);
//  for(i=0; i<num_el; i++){   propagate_electronic(dt_el, Cact, ham.children, prms.rep_tdse, prms.isNBRA);  }
  for(i=0; i<num_el; i++){   
        propagate_electronic(dt_el, Cact, dyn_var.proj_adi, ham, ham_aux, prms, rnd);
  }

  //update_forces(prms, dyn_var, ham);

  //cout<<"Before: "<<dyn_var.act_states[0]<<endl;
  //update_active_states(dyn_var.act_states, dyn_var.proj_adi);
  //cout<<"After: "<<dyn_var.act_states[0]<<endl;

  dyn_var.update_density_matrix(prms, ham, 1); // This one is okay
  dyn_var.act_states = update_active_states(dyn_var.act_states, dyn_var.proj_adi); 
  update_forces(prms, dyn_var, ham);


  // NVT dynamics
  if(prms.ensemble==1){  
    for(traj=0; traj<ntraj; traj++){
      t2[0] = traj; 
      pop_submatrix(p, p_traj, t1, t2);
      double ekin = compute_kinetic_energy(p_traj, invM, prms.thermostat_dofs);
      therm[traj].propagate_nhc(prms.dt, ekin, 0.0, 0.0);
    }

  }

//  F_eff = aux_get_forces(prms, dyn_var, ham); 
//  p = p + F_eff * 0.5 * prms.dt;
  *dyn_var.p = *dyn_var.p + 0.5*prms.dt* (*dyn_var.f);


  // Kinetic constraint
  for(cdof=0; cdof<prms.constrained_dofs.size(); cdof++){   
    p.scale(prms.constrained_dofs[cdof], -1, 0.0); 
  }

  // NVT dynamics
  if(prms.ensemble==1){  
    for(idof=0; idof<n_therm_dofs; idof++){
      dof = prms.thermostat_dofs[idof];
      for(traj=0; traj<ntraj; traj++){
        p.scale(dof, traj, therm[traj].vel_scale(0.5*prms.dt));
      }// traj
    }// idof 
  }


  update_Hamiltonian_variables(prms, dyn_var, ham, ham_aux, py_funct, params, 1);


  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
// TEMPORARY COMMENT
  //dyn_var.update_amplitudes(prms, ham);  // Don't do this - then we are fine with the diabatic picture

//  dyn_var.update_density_matrix(prms, ham, 1); // This one is okay


  //============== Begin the TSH part ===================

  //================= Update decoherence rates & times ================
  /// Effectively turn off decoherence effects
  if(prms.decoherence_times_type==-1){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = 0.0;   }
  }
  /// Just use the plain times given from the input, usually the mSDM formalism
  else if(prms.decoherence_times_type==0){
    for(traj=0; traj<ntraj1; traj++){   decoherence_rates[traj] = *prms.decoherence_rates;   }
  }
  /// Compute the dephasing rates according the original energy-based formalism
  else if(prms.decoherence_times_type==1){
    Eadi = get_Eadi(ham); 
    Ekin = dyn_var.compute_kinetic_energies();  
    decoherence_rates = edc_rates(Eadi, Ekin, prms.decoherence_C_param, prms.decoherence_eps_param, prms.isNBRA);       
  }

  else if(prms.decoherence_times_type==2){
    decoherence_rates = schwartz_1(prms, *dyn_var.ampl_adi, ham, *prms.schwartz_decoherence_inv_alpha); 
  }

  else if(prms.decoherence_times_type==3){
    decoherence_rates = schwartz_2(prms, ham, *prms.schwartz_decoherence_inv_alpha); 
  }

  ///== Optionally, apply the dephasing-informed correction ==
  if(prms.dephasing_informed==1){
    Eadi = get_Eadi(ham); 
    MATRIX ave_gaps(*prms.ave_gaps);
    dephasing_informed_correction(decoherence_rates, Eadi, ave_gaps, prms.isNBRA);
  }

  //============ Apply decoherence corrections ==================
  // SDM and alike methods - only in the adiabatic rep
  if(prms.decoherence_algo==0 && prms.rep_tdse==1){
    *dyn_var.ampl_adi = sdm(*dyn_var.ampl_adi, prms.dt, dyn_var.act_states, decoherence_rates, prms.sdm_norm_tolerance, prms.isNBRA);
  }
  // BCSH
  else if(prms.decoherence_algo==3){ 
// TEMPORARY COMMENTS - next 2 lines
//    *dyn_var.reversal_events = wp_reversal_events(p, invM, act_states, ham, projectors, prms.dt);
//    Coeff = bcsh(Coeff, prms.dt, act_states, *dyn_var.reversal_events);
  }

  // MFSD
  else if(prms.decoherence_algo==4){
    Cact = mfsd(*dyn_var.p, Cact, *dyn_var.iM, prms.dt, decoherence_rates, ham, rnd, prms.isNBRA);
  }



  //========= Use the resulting amplitudes to do the hopping =======
  //dyn_var.update_amplitudes(prms, ham);
  //dyn_var.update_density_matrix(prms, ham, 1);


  // Adiabatic dynamics
  if(prms.tsh_method==-1){ ;; } 

  // FSSH, GFSH, MSSH
  else if(prms.tsh_method == 0 || prms.tsh_method == 1 || prms.tsh_method == 2){

    // Compute hopping probabilities
    // vector<MATRIX> g( hop_proposal_probabilities(prms, q, p, invM, Coeff, ham, prev_ham_dia) );

    /// Compute hop proposal probabilities from the active state of each trajectory to all other states 
    /// of that trajectory
    vector< vector<double> > g;
    g = hop_proposal_probabilities(prms, dyn_var, ham, prev_ham_dia);


    // Propose new discrete states for all trajectories
    vector<int> prop_states( propose_hops(g, dyn_var.act_states, rnd) );

    
    // Decide if to accept the transitions (and then which)
    // Here, it is okay to use the local copies of the q, p, etc. variables, since we don't change the actual variables
    vector<int> old_states(ntraj); //act_states);
    for(traj=0; traj<ntraj; traj++){ old_states[traj] = act_states[traj]; }
    act_states = accept_hops(prms, q, p, invM, Cact, ham, prop_states, act_states, rnd);    


    // Velocity rescaling: however here we may be changing velocities
    //
    handle_hops_nuclear(prms, q, p, invM, Cact, ham, act_states, old_states);


    if(prms.decoherence_algo==1){
      // Instantaneous decoherence 
      instantaneous_decoherence(Cact, act_states, prop_states, old_states, prms.instantaneous_decoherence_variant, prms.collapse_option);
    } 
    else if(prms.decoherence_algo==2){
      /// Temporarily commented AVA 11/7/2022
      ///apply_afssh(dyn_var, Coeff, act_states, invM, ham, dyn_params, rnd);
    }// AFSSH
        
  }// tsh_method == 0, 1, 2, 4
  
/*
  // DISH
  else if(prms.tsh_method == 3 ){

    /// Advance coherence times
    coherence_time.add(-1, -1, prms.dt);

    /// New version, as of 8/3/2020
    vector<int> old_states(act_states);
    act_states = dish(prms, q, p, invM, Coeff, ham, act_states, coherence_time, decoherence_rates, rnd);

    /// Velocity rescaling
    handle_hops_nuclear(prms, q, p, invM, Coeff, ham, act_states, old_states);


  }// tsh_method == 3
*/
  else{   cout<<"tsh_method == "<<prms.tsh_method<<" is undefined.\nExiting...\n"; exit(0);  }



  project_out_states.clear();


  if(prms.rep_tdse==0){ *dyn_var.ampl_dia = Cact; } //*dyn_var.ampl_dia; }
  else if(prms.rep_tdse==1){ *dyn_var.ampl_adi = Cact; } //*dyn_var.ampl_adi; }



}







}// namespace libdyn
}// liblibra

