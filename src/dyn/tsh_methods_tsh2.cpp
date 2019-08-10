/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_methods_tsh2.cpp
  \brief The file implements a new, experimental SH procedure
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{



void tsh2(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, vector<int>& act_states,
          nHamiltonian& ham, bp::object py_funct, bp::object params, boost::python::dict params1, Random& rnd){

/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] act_states - vector of ntraj indices of the physical states in which each of the trajectories
  initially is (active states). 
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) 
  - its internal variables (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] params1 The Python dictionary containing the control parameters passed to this function
  \param[in] rnd The Random number generator object

  Return: propagates q, p and updates state variables

*/


  /**
    Setup the default values of the control parameters:
  */

  int do_reverse = 1;         ///< 0 - do not revert momenta at the frustrated hops, 1 - do revert the momenta
  int do_rescale = 1;         ///< 0 - do not rescale momenta at the accepted hops, 1 - do rescale the momenta
  int vel_rescale_opt = 0;    ///< How to rescale momenta if the hops are successful:
                              ///<   0 - rescale along the directions of derivative couplings

  int ham_rep = 0; ///< default -- assume the Hamiltonian is first computed in the diabatic representation
                   ///< and then will be transformed to the adiabatic in this function. 

  /**
    Extract the parameters from the input dictionary
  */

  std::string key;
  for(int i=0;i<len(params1.values());i++){
    key = extract<std::string>(params1.keys()[i]);

         if(key=="ham_rep") { ham_rep = extract<int>(params1.values()[i]);   }
    else if(key=="do_reverse") { do_reverse = extract<int>(params1.values()[i]);  }
    else if(key=="do_rescale") { do_rescale = extract<int>(params1.values()[i]);  }

  }


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int traj, dof, i, j;


  CMATRIX** Uprev;   Uprev = new CMATRIX*[ntraj];
  CMATRIX X(1, ham.nadi);


  // ====== Get the previous eigenvectors =====  
  for(traj=0; traj<ntraj; traj++){
    Uprev[traj] = new CMATRIX(ham.nadi, ham.nadi);
    *Uprev[traj] = ham.children[traj]->get_basis_transform();  
  }


  // ====== Evolve momenta =====
  p = p + ham.forces_adi(act_states).real() * 0.5*dt;

  // ====== Evolve positions =====
  // For efficiency, switch to the element-wise multiplication
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 
    }
  }

  // ====== Update the potential =====
  if(ham_rep==0){
    ham.compute_diabatic(py_funct, bp::object(q), params, 1);
    ham.compute_adiabatic(1, 1);
  }
  else if(ham_rep==1){
    ham.compute_adiabatic(py_funct, bp::object(q), params, 1);
  }

  // ====== Evolve momenta =====
  p = p + ham.forces_adi(act_states).real() * 0.5 * dt;


  //======= Compute the surface hopping probabilities in
  // a new way: using the evolved eigenvectors ========

  MATRIX g(ham.nadi,ham.nadi);
  vector<int> fstates(ntraj, 0);  // proposed final states

  for(traj=0; traj<ntraj; traj++){
    
    // X = <psi_i(t)|psi(t+dt)>, where i is the active state

    X = Uprev[traj]->col(act_states[traj]).H() * ham.children[traj]->get_basis_transform();
    for(i=0; i<ham.nadi; i++){  g.set(act_states[traj], i, std::norm(X.get(i))); }

    /// Attempt to hop
    double ksi = rnd.uniform(0.0,1.0);  /// generate random number 
    fstates[traj] = hop(act_states[traj], g, ksi); /// Proposed hop

  }// for trajectories

  for(traj=0; traj<ntraj; traj++){  delete Uprev[traj]; }
  delete Uprev;


  // Hop acceptance/rejection - velocity rescaling
  act_states = apply_transition1(p, invM, ham, act_states, fstates, vel_rescale_opt, do_reverse, do_rescale); // non-physical states


}




void tsh2a(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
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
  int rep_sh = 1;             ///< The representation to run the SH : 0 - diabatic, 1 - adiabatic
  int ham_rep = 0;            ///< 0 - diabatic, 1 - adiabatic
  int tsh_method = 0;         ///< Formula for computing SH probabilities: 0 - FSSH, 1 - GFSH, 2 - MSSH
  int use_boltz_factor = 0;   ///< Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
  double Temperature = 300.0; ///< Temperature of the system
  int do_reverse = 1;         ///< 0 - do not revert momenta at the frustrated hops, 1 - do revert the momenta
  int do_rescale = 1;
  int vel_rescale_opt = 0;    ///< How to rescale momenta if the hops are successful:
                              ///<   0 - rescale along the directions of derivative couplings
                              ///<   1 - rescale in the diabatic basis - don't care about the velocity directions, just a uniform rescaling,
                              ///<   2 - do not rescale, as in the NBRA.


  /**
    Extract the parameters from the input dictionary
  */

  std::string key;
  for(int i=0;i<len(params1.values());i++){
    key = extract<std::string>(params1.keys()[i]);

    if(key=="rep") { rep = extract<int>(params1.values()[i]); }
    else if(key=="rep_sh") { rep_sh = extract<int>(params1.values()[i]);  }
    else if(key=="ham_rep") { ham_rep = extract<int>(params1.values()[i]);   }
    else if(key=="tsh_method") { tsh_method = extract<int>(params1.values()[i]);  }
    else if(key=="use_boltz_factor") { use_boltz_factor = extract<int>(params1.values()[i]);  }
    else if(key=="Temperature") { Temperature = extract<double>(params1.values()[i]);  }
    else if(key=="do_reverse") { do_reverse = extract<int>(params1.values()[i]);  }
    else if(key=="do_rescale") { do_rescale = extract<int>(params1.values()[i]);  }
    else if(key=="vel_rescale_opt") { vel_rescale_opt = extract<int>(params1.values()[i]);  }
  }


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;

  vector<int> perm_t; 
  CMATRIX states(nst, ntraj); // CMATRIX version of "act_states"
  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);
  complex<double> one(1.0, 0.0);

  MATRIX g(nst,nst); /// the matrix of hopping probability
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);

  vector<int> fstates(ntraj,0); 


  CMATRIX** Uprev;   Uprev = new CMATRIX*[ntraj];
  CMATRIX X(ham.nadi, ham.nadi);

  // ====== Get the previous eigenvectors =====  
  for(traj=0; traj<ntraj; traj++){
    Uprev[traj] = new CMATRIX(ham.nadi, ham.nadi);
    *Uprev[traj] = ham.children[traj]->get_basis_transform();  
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
  /// tsh_indx2vec(ham, states, act_states);
  if(rep==0){    p = p + ham.forces_adi(act_states).real() * 0.5 * dt;  }  /// !!! TODO
  else if(rep==1){    p = p + ham.forces_adi(act_states).real() * 0.5 * dt;  }


  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 
    }
  }


  // ====== Update the potential =====
  if(ham_rep==0){
    ham.compute_diabatic(py_funct, bp::object(q), params, 1);
    ham.compute_adiabatic(1, 1);
  }
  else if(ham_rep==1){
    ham.compute_adiabatic(py_funct, bp::object(q), params, 1);
  }

  // Update the Ehrenfest forces for all trajectories
  ///  tsh_indx2vec(ham, states, act_states);
  if(rep==0){    p = p + ham.forces_adi(act_states).real() * 0.5 * dt;  }  /// !!! TODO
  else if(rep==1){    p = p + ham.forces_adi(act_states).real() * 0.5 * dt;  }


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
  ///tsh_physical2internal(ham, istates, act_states);


  for(traj=0; traj<ntraj; traj++){

    X = Uprev[traj]->H() * ham.children[traj]->get_basis_transform();

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
    else{
      cout<<"Error in tsh0: tsh_method can be 0, 1, or 2. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    for(i=0;i<nst;i++){
      g.scale(act_states[traj], i, std::norm(X.get(act_states[traj], i)) ); 
    }

    /// Attempt to hop
    double ksi = rnd.uniform(0.0,1.0);  /// generate random number 
    fstates[traj] = hop(act_states[traj], g, ksi); /// Proposed hop


  }// for traj


  for(traj=0; traj<ntraj; traj++){  delete Uprev[traj]; }
  delete Uprev;


  // Hop acceptance/rejection - velocity rescaling
  act_states = apply_transition1(p, invM, ham, act_states, fstates, vel_rescale_opt, do_reverse, 1); // non-physical states

  // Convert from the internal indexing to the physical
//  tsh_internal2physical(ham, istates, act_states);


}





}// namespace libdyn
}// liblibra

