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


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


namespace bp = boost::python;


MATRIX aux_get_forces(dyn_control_params& prms, nHamiltonian& ham, vector<int>& act_states, CMATRIX& amplitudes){
  /**
    Compute the force depending on the method used
  */

  int ndof = ham.nnucl;
  int ntraj = act_states.size();

  MATRIX F(ndof, ntraj);

  if(prms.force_method==0){  

    // Don't compute forces at all - e.g. in NBRA

  }// NBRA

  else if(prms.force_method==1){  

    // TSH or adiabatic (including excited states)
    // state-specific forces

    if(prms.rep_force==0){  
      // Diabatic 
      F = ham.forces_dia(act_states).real();
    }

    else if(prms.rep_force==1){  
      // Adiabatic 
      F = ham.forces_adi(act_states).real();
    }

  }// TSH && adiabatic

  else if(prms.force_method==2){  

    if(prms.rep_force==0){  
      // Diabatic 
      F = ham.Ehrenfest_forces_dia(amplitudes, 1).real();
    }

    else if(prms.rep_force==1){  
      // Adiabatic 
      F = ham.Ehrenfest_forces_adi(amplitudes, 1).real();
    }
  
  }// Ehrenfest

  return F;
}


MATRIX aux_get_forces(bp::dict params, nHamiltonian& ham, vector<int>& act_states, CMATRIX& amplitudes){

  dyn_control_params prms;
  prms.set_parameters(params);

  return aux_get_forces(prms, ham, act_states, amplitudes);
  
}


void aux_get_transforms(CMATRIX** Uprev, nHamiltonian& ham){

  // For adiabatic representation only:
  // Save the previous orbitals info - in case we need to
  // do either phase correction of state tracking

  int ntraj = ham.children.size();

  for(int traj=0; traj<ntraj; traj++){
    *Uprev[traj] = ham.children[traj]->get_basis_transform();  
  }

}



void do_reordering(dyn_control_params& prms, nHamiltonian& ham,
                   vector<int>& act_states, CMATRIX& C, CMATRIX** Uprev, Random& rnd){

  // Reordering, if needed

  int ntraj = ham.children.size();
  int nst = ham.nadi;

  vector<int> perm_t; 
  vector<int> el_stenc_x(nst, 0); for(int i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 

  CMATRIX* St;   St = new CMATRIX(ham.nadi, ham.nadi);


  for(int traj=0; traj<ntraj; traj++){
    *St = (*Uprev[traj]).H() * ham.children[traj]->get_basis_transform();


    if(prms.state_tracking_algo==1){
        perm_t = get_reordering(*St);
    }
    else if(prms.state_tracking_algo==2){
        CMATRIX Eadi = ham.children[traj]->get_ham_adi();
        perm_t = Munkres_Kuhn(*St, Eadi, prms.MK_alpha, prms.MK_verbosity);
    }
    if(prms.state_tracking_algo==3){
        perm_t = get_stochastic_reordering(*St, rnd);
    }

    // Switch the active states
    act_states[traj] = perm_t[ act_states[traj] ];
          
    el_stenc_y[0] = traj;
    CMATRIX x(ham.nadi, 1); 
    x = C.col(traj);
    x.permute_rows(perm_t);
    push_submatrix(C, x, el_stenc_x, el_stenc_y);

  }// for trajectories

  delete St;

}

void do_phase_correction(dyn_control_params& prms, nHamiltonian& ham, 
                         vector<int>& act_states, CMATRIX& C, CMATRIX** Uprev){

  CMATRIX* phases; phases = new CMATRIX(ham.nadi, 1); 

  int ntraj = ham.children.size();
  int nst = ham.nadi;

  vector<int> el_stenc_x(nst, 0); for(int i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 

  for(int traj=0; traj<ntraj; traj++){

    // Phase correction in U, NAC, and Hvib
    *phases = ham.children[traj]->update_phases(*Uprev[traj], 1);

    // Phase correction in Cadi
    el_stenc_y[0] = traj;
    CMATRIX x(ham.nadi, 1); 
    x = C.col(traj);
    phase_correct_ampl(&x, phases);
    push_submatrix(C, x, el_stenc_x, el_stenc_y);

  }// for traj

  delete phases;

}

void update_Hamiltonian_q(dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  //------ Update the internals of the Hamiltonian object --------
  // We call the external function that would do the calculations
  if(prms.rep_tdse==0){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, bp::object(q), model_params, 1);
    }
  }
  if(prms.rep_tdse==1){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, bp::object(q), model_params, 1);
      ham.compute_adiabatic(1, 1);
    }
    else if(prms.rep_ham==1){
      ham.compute_adiabatic(py_funct, bp::object(q), model_params, 1);
    }
  }

}


void update_Hamiltonian_q(bp::dict prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_q(_prms, q, ham, py_funct, model_params);

}


void update_Hamiltonian_p(dyn_control_params& prms, nHamiltonian& ham, 
                          MATRIX& p, MATRIX& invM){

  /**
    Update of the vibronic Hamiltonian in response to changed p
  */

  // Update NACs and Hvib for all trajectories
  if(prms.rep_tdse==0){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_dia(p, invM, 0, 1);
      ham.compute_hvib_dia(1);
    }
  }
  else if(prms.rep_tdse==1){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_adi(p, invM, 0, 1); 
      ham.compute_hvib_adi(1);
    }
  }
}


void update_Hamiltonian_p(bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_p(_prms, ham, p, invM);

}



CMATRIX transform_amplitudes(int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham){
/**
  This function converts the amplitudes from one representation to another

  The reason: we may be solving TD-SE (computing forces) in one representation
  but compute the hopping probabilities in another one.

  This function assumes we already have the basis transformation matrix in ham object 
  computed/updated

*/

  int nst = C.n_rows;    
  int ntraj = C.n_cols;

  CMATRIX Coeff(nst,ntraj); 

  /// Depending on the basis, select which   
  /// C - the basis in which the electron-nuclear propagation is done
  /// Coeff - the basis in which SH is done

  if(rep_in==0){                  // Input in the diabatic basis
    if(rep_out==0){ Coeff = C; }  // Output in the diabatic basis too
    else if(rep_out==1){  ham.ampl_dia2adi(C, Coeff, 0, 1);  } // Output in the adiabatic basis
  }
  else if(rep_in==1){   // Input in the adiabatic basis 
    if(rep_out==0){  ham.ampl_adi2dia(Coeff, C, 0, 1); }  // Output in the diabatic basis 
    else if(rep_out==1){ Coeff = C;  } // Output in the diabatic basis too
  }

  return Coeff;
}



void do_surface_hopping(dyn_control_params& prms,
              MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, vector<MATRIX>& prev_ham_dia, Random& rnd){


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;


  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);


  MATRIX g(nst,nst); /// the matrix of hopping probability
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);
  vector<int> fstates(ntraj,0); 

  //============== Begin the TSH part ===================

  // Proposed hops probabilities
  for(traj=0; traj<ntraj; traj++){

    nucl_stenc_y[0] = traj;
    el_stenc_y[0] = traj;
    full_id[1] = traj;

    pop_submatrix(C, coeff, el_stenc_x, el_stenc_y);

    if(prms.tsh_method == 0){ // FSSH
      g = compute_hopping_probabilities_fssh(coeff, ham.children[traj], prms.rep_sh,
                                   prms.dt, prms.use_boltz_factor, prms.Temperature);
    }
    else if(prms.tsh_method == 1){ // GFSH
      g = compute_hopping_probabilities_gfsh(coeff, ham.children[traj], prms.rep_sh,
                                   prms.dt, prms.use_boltz_factor, prms.Temperature);
    }
    else if(prms.tsh_method == 2){ // MSSH
      g = compute_hopping_probabilities_mssh(coeff);
    }
    else if(prms.tsh_method == 3){ // LZ
      pop_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);
      g = compute_hopping_probabilities_lz(ham.children[traj], prms.rep_lz, p_traj, invM, prev_ham_dia[traj]);
    }

    else{
      cout<<"Error in tsh1: tsh_method can be -1, 0, 1, 2, or 3. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    //============== Compute the potential states after the proposed hops =======================
    double ksi = rnd.uniform(0.0,1.0);             /// generate random number 
    fstates[traj] = hop(act_states[traj], g, ksi); /// Proposed hop

  }// for traj

  //======== Decide whether to accept the proposed hops and if so what to do with the nuclei 
  act_states = apply_transition1(p, invM, ham, act_states, fstates, prms.vel_rescale_opt, prms.do_reverse, 1); 

}


void do_surface_hopping(bp::dict prms,
              MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, vector<MATRIX>& prev_ham_dia, Random& rnd){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  do_surface_hopping(_prms, q, p, invM, C, act_states, ham, prev_ham_dia, rnd);

}


void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, bp::object py_funct, bp::dict params, bp::dict dyn_params, Random& rnd){

/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x ntraj]  or [ndia x ntraj] matrix containing the electronic coordinates
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

  dyn_control_params prms;
  prms.set_parameters(dyn_params);


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof;

  CMATRIX** Uprev; 
  vector<MATRIX> prev_ham_dia(ntraj, MATRIX(nst, nst));


  if(prms.tsh_method == 3){
    for(traj=0; traj<ntraj; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }
  }

  //============ Update the Hamiltonian object =============
  // In case, we may need phase correction & state reordering
  // prepare the temporary files
  if(prms.rep_tdse==1){      
    if(prms.do_phase_correction || prms.state_tracking_algo > 0){

      Uprev = new CMATRIX*[ntraj];
      for(traj=0; traj<ntraj; traj++){  Uprev[traj] = new CMATRIX(nst, nst);  }

      aux_get_transforms(Uprev, ham);
    }
  }// rep == 1



 
  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  propagate_electronic(0.5*prms.dt, C, ham.children, prms.rep_tdse);   

  //============== Nuclear propagation ===================
 
  p = p + aux_get_forces(prms, ham, act_states, C) * 0.5 * prms.dt;

  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * prms.dt ); 
    }
  }


  // Recompute the matrices at the new geometry and apply any
  // necessary fixes 
  update_Hamiltonian_q(prms, q, ham, py_funct, params);



  // Apply phase correction and state reordering as needed
  if(prms.rep_tdse==1){
    if(prms.state_tracking_algo > 0){
      do_reordering(prms, ham, act_states, C, Uprev, rnd);
    }

    if(prms.do_phase_correction){
      do_phase_correction(prms, ham, act_states, C, Uprev);
    }

    if(prms.do_phase_correction || prms.state_tracking_algo >0 ){
      for(traj=0; traj<ntraj; traj++){  delete Uprev[traj]; }  delete Uprev;
    }

  }// rep_tdse == 1


  // Update the Ehrenfest forces for all trajectories
  //tsh_indx2vec(ham, states, act_states);

  p = p + aux_get_forces(prms, ham, act_states, C) * 0.5 * prms.dt;


  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  update_Hamiltonian_p(prms, ham, p, invM);
  propagate_electronic(0.5*prms.dt, C, ham.children, prms.rep_tdse);   



  //============== Begin the TSH part ===================

  // To be able to compute transition probabilities, compute the corresponding amplitudes
  CMATRIX Coeff(nst,ntraj); 
  Coeff = transform_amplitudes(prms.rep_tdse, prms.rep_sh, C, ham);

  // Use them to do the hopping
  if(prms.tsh_method >=0){
    do_surface_hopping(prms, q, p, invM, Coeff, act_states, ham, prev_ham_dia, rnd);
  }

}


}// namespace libdyn
}// liblibra

