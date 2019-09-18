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


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

struct dyn_control_params{

  /**
    The representation to run the Ehrenfest : 0 - diabatic, 1 - adiabatic
  */
  int rep;


  /** 
   The representation of the Hamiltonian update: 0 - diabatic, 1 - adiabatic
   this is the representation in which the computed properties are assumed to be
   for instance, we may have set it to 1, to read the adiabatic energies and couplings,
   to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
   calculations, or with the NBRA
  */
  int ham_rep;


  /**
   The representation to run the SH : 0 - diabatic, 1 - adiabatic
  */
  int rep_sh;


  /** 
    The representation to compute LZ probabilitieis: 0 - diabatic, 1- adiabatic 
  */
  int rep_lz;


  /** 
    Formula for computing SH probabilities: 0 - FSSH, 1 - GFSH, 2 - MSSH
  */
  int tsh_method;


  /** 
    How to compute forces in the dynamics: 
      0 - don't compute forces at all - e.g. we do not really need them
      1 - state-specific  as in the TSH or adiabatic (including adiabatic excited states)
      2 - Ehrenfest
  */
  int force_method;

  /** 
    How to update NACs and vibronic Hamiltonian before electronic TD-SE propagation
      0 - don't update them (e.g. for simplest NAC)
      1 - update according to changed momentum and existing derivative couplings
  */
  int nac_update_method;


  /** 
    In which representation to copute forces: 
      0 - diabatic
      1 - adiabatic
  */
  int rep_force;


  /**
    Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale                            
  */
  int use_boltz_factor;

  
  /**
    Temperature of the system
  */ 
  double Temperature;


  /** 
    Do not revert momenta at the frustrated hops, 1 - do revert the momenta
  */
  int do_reverse;


  /** 
    How to rescale momenta if the hops are successful:
      0 - rescale along the directions of derivative couplings
      1 - rescale in the diabatic basis - don't care about the velocity directions, just a uniform rescaling,
      2 - do not rescale, as in the NBRA.  
  */
  int vel_rescale_opt;


  /** 
    integration timestep [units: a.u., default: 41 a.u. = 1 fs]
  */
  double dt;

  /** 
    Option to perform the phase correction: 0 - no, 1 - yes (default)
  */
  int do_phase_correction;


  /** 
    State tracking algorithm:
      0 - no state tracking
      1 - method of Kosuke Sato (may fail by getting trapped into an infinite loop)
      2 - Munkres-Kuhn (Hungarian) algorithm (default)
  */
  int state_tracking_algo;

  /** 
    Munkres-Kuhn alpha (selects the range of orbitals included in reordering) [default: 0.0]
  */
  double MK_alpha;

  /**
    Munkres-Kuhn verbosity: 0 - no extra output (default), 1 - details
  */
  int MK_verbosity;


  /**
    Active electronic state used in the adiabatic MD: 0 - ground state, 1 - first excited, and so on
  */
  int act_state; 


  /**
    A selector of a method to couple the trajectories in this ensemble:
      0 - no coupling, 1 - ETHD, 2 - ETHD3 (experimental), 22 - another flavor of ETHD3 (experimental)
  */
  int entanglement_opt;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 )
  */
  double ETHD3_alpha;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the momentum space, that is   ~exp(-beta*(P-P0)^2 )
  */
  double ETHD3_beta;



};

void set_defaults(dyn_control_params& prms){
/**

  This function initializes the default values of control parameters

*/

  prms.rep = 1;
  prms.ham_rep = 0;
  prms.rep_sh = 1;
  prms.rep_lz = 0;
  prms.tsh_method = 0;
  prms.force_method = 1;
  prms.nac_update_method = 1;
  prms.rep_force = 1;
  prms.use_boltz_factor = 0;
  prms.Temperature = 300.0;
  prms.do_reverse = 1;
  prms.vel_rescale_opt = 0;
  prms.dt = 41.0;
  prms.do_phase_correction = 1;
  prms.state_tracking_algo = 2;
  prms.MK_alpha = 0.0;
  prms.MK_verbosity = 0;

  prms.act_state = 0; 

  prms.entanglement_opt = 0;
  prms.ETHD3_alpha = 1.0;
  prms.ETHD3_beta = 1.0;


}

void set_parameters(dyn_control_params& prms, boost::python::dict params){
/**
  Extract the parameters from the input dictionary
*/

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = extract<std::string>(params.keys()[i]);


    if(key=="rep") { prms.rep = extract<int>(params.values()[i]); }
    else if(key=="ham_rep") { prms.ham_rep = extract<int>(params.values()[i]);   }
    else if(key=="rep_sh") { prms.rep_sh = extract<int>(params.values()[i]);  }
    else if(key=="rep_lz") { prms.rep_lz = extract<int>(params.values()[i]);  }
    else if(key=="tsh_method") { prms.tsh_method = extract<int>(params.values()[i]);  }
    else if(key=="force_method") { prms.force_method = extract<int>(params.values()[i]);  }
    else if(key=="nac_update_method") { prms.nac_update_method = extract<int>(params.values()[i]);  }
    else if(key=="rep_force") { prms.rep_force = extract<int>(params.values()[i]);  }
    else if(key=="use_boltz_factor") { prms.use_boltz_factor = extract<int>(params.values()[i]);  }
    else if(key=="Temperature") { prms.Temperature = extract<double>(params.values()[i]);  }
    else if(key=="do_reverse") { prms.do_reverse = extract<int>(params.values()[i]);  }
    else if(key=="vel_rescale_opt") { prms.vel_rescale_opt = extract<int>(params.values()[i]);  }

    else if(key=="dt") { prms.dt = extract<double>(params.values()[i]);  }

    // Phase correction
    else if(key=="do_phase_correction") { prms.do_phase_correction = extract<int>(params.values()[i]);  }

    // State tracking options
    else if(key=="state_tracking_algo"){  prms.state_tracking_algo = extract<int>(params.values()[i]);  }
    else if(key=="MK_alpha") { prms.MK_alpha = extract<double>(params.values()[i]);  }
    else if(key=="MK_verbosity") { prms.MK_verbosity = extract<int>(params.values()[i]);  }

    // Adiabatic dynamics 
    else if(key=="act_state"){ prms.act_state = extract<int>(params.values()[i]); }

    // Trajectory coupling
    else if(key=="entanglement_opt"){ prms.entanglement_opt = extract<int>(params.values()[i]); }
    else if(key=="ETHD3_alpha") { prms.ETHD3_alpha = extract<double>(params.values()[i]);   }
    else if(key=="ETHD3_beta") { prms.ETHD3_beta = extract<double>(params.values()[i]);   }


  }

}

void sanity_check(dyn_control_params& prms){

  if(prms.state_tracking_algo==0 || prms.state_tracking_algo==1 ||
     prms.state_tracking_algo==2 || prms.state_tracking_algo==3){ ; ; }
  else{
    cout<<"Error in tsh1: state_tracking_algo = "<<prms.state_tracking_algo<<" is not allowed\n";
    exit(0);
  }

}


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



void allocate_tmp(CMATRIX** Uprev, nHamiltonian& ham){

  // For adiabatic representation only:
  // Save the previous orbitals info - in case we need to
  // do either phase correction of state tracking

  int ntraj = ham.children.size();

  Uprev = new CMATRIX*[ntraj];
  for(int traj=0; traj<ntraj; traj++){
    Uprev[traj] = new CMATRIX(ham.nadi, ham.nadi);
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
                          vector<int>& act_states, CMATRIX& C, 
                          bp::object py_funct, bp::object params, Random& rnd){

  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  CMATRIX** Uprev; 
  int ntraj = C.n_cols;    
  int traj;


  //============ Update the Hamiltonian object =============
  // In case, we may need phase correction & state reordering
  // prepare the temporary files
  if(prms.rep==1){      
    if(prms.do_phase_correction || prms.state_tracking_algo >0){
      allocate_tmp(Uprev, ham);
    }
  }// rep == 1


  //------ Update the internals of the Hamiltonian object --------
  // We call the external function that would do the calculations
  if(prms.rep==0){      
    if(prms.ham_rep==0){
      ham.compute_diabatic(py_funct, bp::object(q), params, 1);
    }
  }
  if(prms.rep==1){      
    if(prms.ham_rep==0){
      ham.compute_diabatic(py_funct, bp::object(q), params, 1);
      ham.compute_adiabatic(1, 1);
    }
    else if(prms.ham_rep==1){
      ham.compute_adiabatic(py_funct, bp::object(q), params, 1);
    }
  }


  // Apply phase correction and state reordering as needed
  if(prms.rep==1){
    if(prms.state_tracking_algo > 0){
      do_reordering(prms, ham, act_states, C, Uprev, rnd);
    }

    if(prms.do_phase_correction){
      do_phase_correction(prms, ham, act_states, C, Uprev);
    }

    if(prms.do_phase_correction || prms.state_tracking_algo >0 ){
      for(traj=0; traj<ntraj; traj++){  delete Uprev[traj]; }
      delete Uprev;
    }

  }// rep == 1

}


void update_Hamiltonian_p(dyn_control_params& prms, nHamiltonian& ham, 
                          MATRIX& p, MATRIX& invM){

  /**
    Update of the vibronic Hamiltonian in response to changed p
  */

  // Update NACs and Hvib for all trajectories
  if(prms.rep==0){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_dia(p, invM, 0, 1);
      ham.compute_hvib_dia(1);
    }
  }
  else if(prms.rep==1){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_adi(p, invM, 0, 1); 
      ham.compute_hvib_adi(1);
    }

  }

}


CMATRIX transform_amplitudes(dyn_control_params& prms, CMATRIX& C, nHamiltonian& ham){
/**
  This function converts the amplitudes from one representation to another

  The reason: we may be solving TD-SE (computing forces) in one representation
  but compute the hopping probabilities in another one.

*/

  int nst = C.n_rows;    
  int ntraj = C.n_cols;

  CMATRIX Coeff(nst,ntraj); 

  /// Depending on the basis, select which   
  /// C - the basis in which the electron-nuclear propagation is done
  /// Coeff - the basis in which SH is done

  if(prms.rep==0){   // Propagation in the diabatic basis (C - is diabatic)
    if(prms.rep_sh==0){ Coeff = C; }  // SH in the diabatic basis (Coeff - diabatic)
    else if(prms.rep_sh==1){  ham.ampl_dia2adi(C, Coeff, 0, 1);  } // SH in the adiabatic basis (Coeff - adiabatic)
  }
  else if(prms.rep==1){   // Propagation in the adiabatic basis (C - is adiabatic)
    if(prms.rep_sh==0){  ham.ampl_adi2dia(Coeff, C, 0, 1); }  // SH in the diabatic basis (Coeff - diabatic)
    else if(prms.rep_sh==1){ Coeff = C;  } // SH in the adiabatic basis (Coeff - adiabatic)
  }

  return Coeff;
}


void do_surface_hopping(dyn_control_params& prms,
              MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, bp::object py_funct, bp::object params, boost::python::dict params1, 
              vector<MATRIX>& prev_ham_dia, Random& rnd){


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
      cout<<"Error in tsh1: tsh_method can be 0, 1, 2, or 3. Other values are not defined\n";
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

void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, bp::object py_funct, bp::object params, boost::python::dict params1, Random& rnd){

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

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof;


  dyn_control_params prms;
  set_defaults(prms);
  set_parameters(prms, params1);
  sanity_check(prms);


  vector<MATRIX> prev_ham_dia(ntraj, MATRIX(nst, nst));

  if(prms.tsh_method == 3){
    for(traj=0; traj<ntraj; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }
  }

 
  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  propagate_electronic(0.5*prms.dt, C, ham.children, prms.rep);   

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
  update_Hamiltonian_q(prms, q, ham, act_states, C, py_funct, params, rnd);

  // Update the Ehrenfest forces for all trajectories
  //tsh_indx2vec(ham, states, act_states);

  p = p + aux_get_forces(prms, ham, act_states, C) * 0.5 * prms.dt;


  //============== Electronic propagation ===================
  // Evolve electronic DOFs for all trajectories
  update_Hamiltonian_p(prms, ham, p, invM);
  propagate_electronic(0.5*prms.dt, C, ham.children, prms.rep);   



  //============== Begin the TSH part ===================

  // To be able to compute transition probabilities, compute the corresponding amplitudes
  CMATRIX Coeff(nst,ntraj); 
  Coeff = transform_amplitudes(prms, C, ham);

  // Use them to do the hopping
  do_surface_hopping(prms, q, p, invM, Coeff, act_states, ham, py_funct, params, params1, prev_ham_dia, rnd);


}


}// namespace libdyn
}// liblibra

