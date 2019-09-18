/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_methods_tsh.cpp
  \brief The file implements the basic (Tully's) surface hopping approach.
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{



void tsh1b(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
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


  /**
    Setup the default values of the control parameters:
  */

  int rep = 1;                ///< The representation to run the Ehrenfest : 0 - diabatic, 1 - adiabatic
  int ham_rep = 0;            ///< The representation of the Hamiltonian update: 0 - diabatic, 1 - adiabatic
                              ///< this is the representation in which the computed properties are assumed to be
                              ///< For instance, we may have set it to 1, to read the adiabatic energies and couplings,
                              ///< to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
                              ///< calculations, or with the NBRA
  int rep_sh = 1;             ///< The representation to run the SH : 0 - diabatic, 1 - adiabatic
  int rep_lz = 0;             ///< The representation to compute LZ probabilitieis: 0 - diabatic, 1- adiabatic 
  int tsh_method = 0;         ///< Formula for computing SH probabilities: 0 - FSSH, 1 - GFSH, 2 - MSSH
  int use_boltz_factor = 0;   ///< Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
  double Temperature = 300.0; ///< Temperature of the system
  int do_reverse = 1;         ///< 0 - do not revert momenta at the frustrated hops, 1 - do revert the momenta
  int vel_rescale_opt = 0;    ///< How to rescale momenta if the hops are successful:
                              ///<   0 - rescale along the directions of derivative couplings
                              ///<   1 - rescale in the diabatic basis - don't care about the velocity directions, just a uniform rescaling,
                              ///<   2 - do not rescale, as in the NBRA.

   double dt = 41.0;          ///< integration timestep [units: a.u., default: 41 a.u. = 1 fs]

   int do_phase_correction = 1; /// Option to perform the phase correction: 0 - no, 1 - yes (default)
   int state_tracking_algo = 2;   ///< State tracking algorithm:
                              ///<   0 - no state tracking
                              ///<   1 - method of Kosuke Sato (may fail by getting trapped into an infinite loop)
                              ///<   2 - Munkres-Kuhn (Hungarian) algorithm (default)
   double MK_alpha = 0.0;     ///< Munkres-Kuhn alpha (selects the range of orbitals included in reordering) [default: 0.0]
   int MK_verbosity = 0;      ///< Munkres-Kuhn verbosity: 0 - no extra output (default), 1 - details


  /**
    Extract the parameters from the input dictionary
  */


  std::string key;
  for(int i=0;i<len(params1.values());i++){
    key = extract<std::string>(params1.keys()[i]);

    if(key=="rep") { rep = extract<int>(params1.values()[i]); }
    else if(key=="ham_rep") { ham_rep = extract<int>(params1.values()[i]);   }
    else if(key=="rep_sh") { rep_sh = extract<int>(params1.values()[i]);  }
    else if(key=="rep_lz") { rep_lz = extract<int>(params1.values()[i]);  }
    else if(key=="tsh_method") { tsh_method = extract<int>(params1.values()[i]);  }
    else if(key=="use_boltz_factor") { use_boltz_factor = extract<int>(params1.values()[i]);  }
    else if(key=="Temperature") { Temperature = extract<double>(params1.values()[i]);  }
    else if(key=="do_reverse") { do_reverse = extract<int>(params1.values()[i]);  }
    else if(key=="vel_rescale_opt") { vel_rescale_opt = extract<int>(params1.values()[i]);  }


    else if(key=="dt") { dt = extract<double>(params1.values()[i]);  }

    // Phase correction
    else if(key=="do_phase_correction") { do_phase_correction = extract<int>(params1.values()[i]);  }

    // State tracking options
    else if(key=="state_tracking_algo"){  state_tracking_algo = extract<int>(params1.values()[i]);  }
    else if(key=="MK_alpha") { MK_alpha = extract<double>(params1.values()[i]);  }
    else if(key=="MK_verbosity") { MK_verbosity = extract<int>(params1.values()[i]);  }

  }

  if(state_tracking_algo==0 || state_tracking_algo==1 || state_tracking_algo==2 || state_tracking_algo==3){ ; ; }
  else{
    cout<<"Error in tsh1: state_tracking_algo = "<<state_tracking_algo<<" is not allowed\n";
    exit(0);
  }


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;

  CMATRIX** Uprev; 
  CMATRIX* X;
  vector<int> perm_t; 
  CMATRIX states(nst, ntraj); // CMATRIX version of "act_states"
  vector<MATRIX> prev_ham_dia(ntraj, MATRIX(nst, nst));

  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);
  complex<double> one(1.0, 0.0);

  MATRIX g(nst,nst); /// the matrix of hopping probability
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);

  vector<int> istates(ntraj,0); 
  vector<int> fstates(ntraj,0); 


  if(tsh_method == 3){
    for(traj=0; traj<ntraj; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }
  }
 
  //============== Electronic propagation ===================
  // Update NACs and Hvib for all trajectories
  if(rep==0){  
    ham.compute_nac_dia(p, invM, 0, 1);
    ham.compute_hvib_dia(1);
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM, 0, 1); 
    ham.compute_hvib_adi(1);
  }

  // Evolve electronic DOFs for all trajectories
  propagate_electronic(0.5*dt, C, ham.children, rep);   

  //============== Nuclear propagation ===================

  // Update the Ehrenfest forces for all trajectories
  tsh_indx2vec(ham, states, act_states);
  if(rep==0){  p = p + ham.Ehrenfest_forces_dia(states, 1).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(states, 1).real() * 0.5*dt;  }


  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 
    }
  }


  // For adiabatic representation only:
  // Save the previous orbitals info - in case we need to
  // do either phase correction of state tracking
  if(rep==1){      
    if(do_phase_correction || state_tracking_algo >0){
      Uprev = new CMATRIX*[ntraj];
      for(traj=0; traj<ntraj; traj++){
        Uprev[traj] = new CMATRIX(ham.nadi, ham.nadi);
        *Uprev[traj] = ham.children[traj]->get_basis_transform();  
      }
    }
  }// rep == 1


  /*
  ham.compute_diabatic(py_funct, bp::object(q), params, 1);
  ham.compute_adiabatic(1, 1);
  */

  if(ham_rep==0){
    ham.compute_diabatic(py_funct, bp::object(q), params, 1);
    ham.compute_adiabatic(1, 1);
  }
  else if(ham_rep==1){
    ham.compute_adiabatic(py_funct, bp::object(q), params, 1);
  }



  if(rep==1){

    // Reordering, if needed
    if(state_tracking_algo > 0){

      X = new CMATRIX(ham.nadi, ham.nadi);

      for(traj=0; traj<ntraj; traj++){
        *X = (*Uprev[traj]).H() * ham.children[traj]->get_basis_transform();


        if(state_tracking_algo==1){
            perm_t = get_reordering(*X);
        }
        else if(state_tracking_algo==2){
            CMATRIX Eadi = ham.children[traj]->get_ham_adi();
            perm_t = Munkres_Kuhn(*X, Eadi, MK_alpha, MK_verbosity);
        }
        if(state_tracking_algo==3){
            perm_t = get_stochastic_reordering(*X, rnd);
        }


        //ham.children[traj]->update_ordering(perm_t, 1);

        
        // Switch the active states
        act_states[traj] = perm_t[ act_states[traj] ];
        
        
        el_stenc_y[0] = traj;
        CMATRIX x(ham.nadi, 1); 
        x = C.col(traj);
        x.permute_rows(perm_t);
        push_submatrix(C, x, el_stenc_x, el_stenc_y);



      }// for trajectories

      delete X;
    }// do_reordering


    if(do_phase_correction){

      CMATRIX* phases; phases = new CMATRIX(ham.nadi, 1); 

      for(traj=0; traj<ntraj; traj++){

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
    }// phase correction

    if(do_phase_correction || state_tracking_algo >0 ){
      for(traj=0; traj<ntraj; traj++){  delete Uprev[traj]; }
      delete Uprev;
    }

  }// rep == 1


  // Update the Ehrenfest forces for all trajectories
  //tsh_indx2vec(ham, states, act_states);

  if(rep==0){  p = p + ham.Ehrenfest_forces_dia(states, 1).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(states, 1).real() * 0.5*dt;  }  



  //============== Electronic propagation ===================
  // Update NACs and Hvib for all trajectories
  if(rep==0){  
    ham.compute_nac_dia(p, invM, 0, 1);
    ham.compute_hvib_dia(1);
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM, 0, 1); 
    ham.compute_hvib_adi(1);
  }

  // Evolve electronic DOFs for all trajectories
  propagate_electronic(0.5*dt, C, ham.children, rep);   

  //============== Begin the TSH part ===================

  CMATRIX Coeff(nst,ntraj); 

  /// Depending on the basis, select which   
  /// C - the basis in which the electron-nuclear propagation is done
  /// Coeff - the basis in which SH is done

  if(rep==0){   // Propagation in the diabatic basis (C - is diabatic)
    if(rep_sh==0){  Coeff = C; }  // SH in the diabatic basis (Coeff - diabatic)
    else if(rep_sh==1){ ham.ampl_dia2adi(C, Coeff, 0, 1);  } // SH in the adiabatic basis (Coeff - adiabatic)
  }
  else if(rep==1){   // Propagation in the adiabatic basis (C - is adiabatic)
    if(rep_sh==0){  ham.ampl_adi2dia(Coeff, C, full_id); }  // SH in the diabatic basis (Coeff - diabatic)
    else if(rep_sh==1){ Coeff = C;  } // SH in the adiabatic basis (Coeff - adiabatic)
  }


  /// Compute the proposed multi-trajectory states
  //istates = tsh_vec2indx(states);  /// starting (non-phisical!) states
//  tsh_physical2internal(ham, istates, act_states);


  for(traj=0; traj<ntraj; traj++){
    nucl_stenc_y[0] = traj;
    el_stenc_y[0] = traj;
    full_id[1] = traj;

    pop_submatrix(Coeff, coeff, el_stenc_x, el_stenc_y);

    if(tsh_method == 0){ // FSSH
      g = compute_hopping_probabilities_fssh(coeff, ham.children[traj], rep_sh, dt, use_boltz_factor, Temperature);
    }
    else if(tsh_method == 1){ // GFSH
      g = compute_hopping_probabilities_gfsh(coeff, ham.children[traj], rep_sh, dt, use_boltz_factor, Temperature);
    }
    else if(tsh_method == 2){ // MSSH
      g = compute_hopping_probabilities_mssh(coeff);
    }
    else if(tsh_method == 3){ // LZ
      pop_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);
      g = compute_hopping_probabilities_lz(ham.children[traj], rep_lz, p_traj, invM, prev_ham_dia[traj]);
    }

    else{
      cout<<"Error in tsh1: tsh_method can be 0, 1, 2, or 3. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    /// Attempt to hop
    double ksi = rnd.uniform(0.0,1.0);  /// generate random number 
    fstates[traj] = hop(act_states[traj], g, ksi); /// Proposed hop
  }// for traj


  // Hop acceptance/rejection - velocity rescaling
  act_states = apply_transition1(p, invM, ham, act_states, fstates, vel_rescale_opt, do_reverse, 1); // non-physical states

  // Convert from the internal indexing to the physical
  // tsh_internal2physical(ham, istates, act_states);


/*
  // Update the Ehrenfest forces for all trajectories
  tsh_indx2vec(ham, states, act_states);
  if(rep==0){  p = p + ham.Ehrenfest_forces_dia(states, 1).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(states, 1).real() * 0.5*dt;  }
*/



}


}// namespace libdyn
}// liblibra

