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
  \param[in,out] C nadi x nadi or ndia x ndia matrix containing the electronic coordinates
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

  int rep = 0;                ///< The representation to run the Ehrenfest : 0 - diabatic, 1 - adiabatic
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
  CMATRIX D(nst,nst);

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



}




}// namespace libdyn
}// liblibra

