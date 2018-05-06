/*********************************************************************************
* Copyright (C) 2015-2018 Alexey V. Akimov
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


int tsh0(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, int state,
         nHamiltonian& ham, bp::object py_funct, bp::object params,  boost::python::dict params1, Random& rnd){

/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x 1] or [ndia x 1] matrix containing the electronic coordinates
  \param[in] state is the index of the currently occupied state (active state)
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] params1 The Python dictionary containing the control parameters passed to this function
  \param[in] rnd The Random number generator object

  Return: the index of the active state at the end of the propagation.

*/



  /**
    Setup the default values of the control parameters:
  */

  int rep = 1;                ///< The representation to run the Ehrenfest : 0 - diabatic, 1 - adiabatic
  int rep_sh = 1;             ///< The representation to run the SH : 0 - diabatic, 1 - adiabatic
  int tsh_method = 0;         ///< Formula for computing SH probabilities: 0 - FSSH, 1 - GFSH, 2 - MSSH
  int use_boltz_factor = 0;   ///< Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
  double Temperature = 300.0; ///< Temperature of the system
  int do_reverse = 1;         ///< 0 - do not revert momenta at the frustrated hops, 1 - do revert the momenta
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
    else if(key=="tsh_method") { tsh_method = extract<int>(params1.values()[i]);  }
    else if(key=="use_boltz_factor") { use_boltz_factor = extract<int>(params1.values()[i]);  }
    else if(key=="Temperature") { Temperature = extract<double>(params1.values()[i]);  }
    else if(key=="do_reverse") { do_reverse = extract<int>(params1.values()[i]);  }
    else if(key=="vel_rescale_opt") { vel_rescale_opt = extract<int>(params1.values()[i]);  }
  }




  int ndof = q.n_rows;
  int dof;

  int nst = 0;
  if(rep==0){  nst = ham.ndia;  }
  else if(rep==1){  nst = ham.nadi;  }

  CMATRIX cstate(nst, 1);
  cstate.set(state, 0, complex<double>(1.0, 0.0));

 
  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   

  //============== Nuclear propagation ===================
    
       if(rep==0){  p = p + ham.forces_dia(cstate).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.forces_adi(cstate).real() * 0.5*dt;  }


  for(dof=0; dof<ndof; dof++){  
    q.add(dof, 0,  invM.get(dof,0) * p.get(dof,0) * dt ); 
  }

  ham.compute_diabatic(py_funct, bp::object(q), params);
  ham.compute_adiabatic(1);


       if(rep==0){  p = p + ham.forces_dia(cstate).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.forces_adi(cstate).real() * 0.5*dt;  }

  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   


  //============== Begin the TSH part ===================

  MATRIX g(nst,nst); /// the matrix of hopping probability

  /// Depending on the basis, select which 
  CMATRIX D(nst,1);

  if(rep==0){   // Propagation in the diabatic basis
    if(rep_sh==0){  D = C; }  // SH in the diabatic basis
    else if(rep_sh==1){ ham.ampl_dia2adi(C, D);  } // SH in the adiabatic basis
  }
  else if(rep==1){   // Propagation in the adiabatic basis
    if(rep_sh==0){  ham.ampl_adi2dia(D, C); }  // SH in the diabatic basis
    else if(rep_sh==1){ D = C;  } // SH in the adiabatic basis
  }


  /// Compute hopping probabilities
  if(tsh_method == 0){ // FSSH
    g = compute_hopping_probabilities_fssh(D, ham, rep_sh, dt, use_boltz_factor, Temperature);
  }
  else if(tsh_method == 1){ // GFSH
    g = compute_hopping_probabilities_gfsh(D, ham, rep_sh, dt, use_boltz_factor, Temperature);
  }
  else if(tsh_method == 2){ // MSSH
    g = compute_hopping_probabilities_mssh(D);
  }
  else{
    cout<<"Error in tsh0: tsh_method can be 0, 1, or 2. Other values are not defined\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  /// Attempt to hop
  double ksi = rnd.uniform(0.0,1.0);  /// generate random number 
  int new_state = hop(state, g, ksi); /// Proposed hop

  /// Check whether the proposed hop should be accepted.
  /// If this it is: the nuclear momenta will be re-scaled

  if(vel_rescale_opt==0){
    new_state = rescale_velocities_adiabatic(p, invM, ham, new_state, state, do_reverse);
  }
  else if(vel_rescale_opt==1){
    new_state = rescale_velocities_diabatic(p, invM, ham, new_state, state);
  }
  else if(vel_rescale_opt==2){
   ;;  /// Don't do anything extra
  }

  return new_state;


}






void tsh1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, CMATRIX& states,
         nHamiltonian& ham, bp::object py_funct, bp::object params, boost::python::dict params1, Random& rnd){

/**
  \brief One step of the TSH algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] C [nadi x ntraj]  or [ndia x ntraj] matrix containing the electronic coordinates
  \param[in,out] states - [nadi x ntraj] or [ndia x ntraj] matrix containing the indices of the states 
  currently occupied by each of the trajectories (active states). Each column of the matrix corresponds to a
  trajectory and contains only 0 in all elements but one (the active state). In that place, it contains 1.
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
  int tsh_method = 0;         ///< Formula for computing SH probabilities: 0 - FSSH, 1 - GFSH, 2 - MSSH
  int use_boltz_factor = 0;   ///< Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
  double Temperature = 300.0; ///< Temperature of the system
  int do_reverse = 1;         ///< 0 - do not revert momenta at the frustrated hops, 1 - do revert the momenta
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
    else if(key=="tsh_method") { tsh_method = extract<int>(params1.values()[i]);  }
    else if(key=="use_boltz_factor") { use_boltz_factor = extract<int>(params1.values()[i]);  }
    else if(key=="Temperature") { Temperature = extract<double>(params1.values()[i]);  }
    else if(key=="do_reverse") { do_reverse = extract<int>(params1.values()[i]);  }
    else if(key=="vel_rescale_opt") { vel_rescale_opt = extract<int>(params1.values()[i]);  }
  }


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = C.n_rows;    
  int traj, dof, i;


  vector<int> nucl_stenc_x(ndof, 0); for(i=0;i<ndof;i++){  nucl_stenc_x[i] = i; }
  vector<int> nucl_stenc_y(1, 0); 
  vector<int> el_stenc_x(nst, 0); for(i=0;i<nst;i++){  el_stenc_x[i] = i; }
  vector<int> el_stenc_y(1, 0); 
  vector<int> full_id(2,0);
  complex<double> one(1.0, 0.0);

 
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
  if(rep==0){  p = p + ham.Ehrenfest_forces_dia(states, 1).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(states, 1).real() * 0.5*dt;  }


  // Update coordinates of nuclei for all trajectories
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  
      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 
    }
  }

  ham.compute_diabatic(py_funct, bp::object(q), params, 1);
  ham.compute_adiabatic(1, 1);


  // Update the Ehrenfest forces for all trajectories
  if(rep==0){  p = p + ham.Ehrenfest_forces_dia(C, 1).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(C, 1).real() * 0.5*dt;  }


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


//  exit(0);

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

  MATRIX g(nst,nst); /// the matrix of hopping probability
  MATRIX p_traj(ndof, 1);
  CMATRIX coeff(nst, 1);

  /// Compute hopping probabilities and exectue the hops
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
    else{
      cout<<"Error in tsh0: tsh_method can be 0, 1, or 2. Other values are not defined\n";
      cout<<"Exiting...\n";
      exit(0);
    }

    /// Attempt to hop
    double ksi = rnd.uniform(0.0,1.0);  /// generate random number 
    int istate;
    complex<double> tmp;
    states.max_col_elt(traj, tmp, istate);
    int new_state = hop(istate, g, ksi); /// Proposed hop

    /// Check whether the proposed hop should be accepted.
    /// If this it is: the nuclear momenta will be re-scaled

    pop_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);

    if(vel_rescale_opt==0){
      new_state = rescale_velocities_adiabatic(p_traj, invM, ham.children[traj], new_state, istate, do_reverse);
    }
    else if(vel_rescale_opt==1){
      new_state = rescale_velocities_diabatic(p_traj, invM, ham.children[traj], new_state, istate);
    }
    else if(vel_rescale_opt==2){
     ;;  /// Don't do anything extra 
    }

    states.set(istate, traj, complex<double>(0.0, 0.0));
    states.set(new_state, traj, complex<double>(1.0, 0.0));

    push_submatrix(p, p_traj, nucl_stenc_x, nucl_stenc_y);

  }// for traj


}



}// namespace libdyn
}// liblibra

