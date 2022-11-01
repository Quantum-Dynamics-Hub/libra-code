/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_ham.cpp
  \brief The file implements the methods that update internal Hamiltonian variables in response to
  change of the dynamical variables
*/


#include "dyn_ham.h"

/// liblibra namespace
namespace liblibra{

using namespace libcalculators;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;


void update_Hamiltonian_q(dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  //------ Update the internals of the Hamiltonian object --------
  // We call the external function that would do the calculations
  if(prms.rep_tdse==0){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, q, model_params, 1);
    }
  }
  if(prms.rep_tdse==1){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, q, model_params, 1);
      ham.compute_adiabatic(1, 1);
    }
    else if(prms.rep_ham==1){
      ham.compute_adiabatic(py_funct, q, model_params, 1);
    }
  }

}


void update_Hamiltonian_q(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  //------ Update the internals of the Hamiltonian object --------
  // We call the external function that would do the calculations
  if(prms.rep_tdse==0){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, *dyn_var.q, model_params, 1);
    }
  }
  if(prms.rep_tdse==1){      
    if(prms.rep_ham==0){
      ham.compute_diabatic(py_funct, *dyn_var.q, model_params, 1);
      ham.compute_adiabatic(1, 1);
    }
    else if(prms.rep_ham==1){
      ham.compute_adiabatic(py_funct, *dyn_var.q, model_params, 1);
    }
  }

}



void update_Hamiltonian_q(bp::dict prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_q(_prms, q, ham, py_funct, model_params);

}

void update_Hamiltonian_q(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_q(_prms, *dyn_var.q, ham, py_funct, model_params);

}



void update_Hamiltonian_q_ethd(dyn_control_params& prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params, MATRIX& invM){

//  cout<<"updating_Hamiltonian_q_ethd\n Option = "<<prms.entanglement_opt<<endl;

  if(prms.entanglement_opt==0){    /* Nothing to do */   }
  else if(prms.entanglement_opt==1){   ham.add_ethd_adi(q, invM, 1);  }
  else if(prms.entanglement_opt==2){   ham.add_ethd3_adi(q, invM, prms.ETHD3_alpha, 1);  }
  else if(prms.entanglement_opt==22){  ham.add_ethd3_adi(q, p, invM, prms.ETHD3_alpha, prms.ETHD3_beta, 1);  }
  else{
    cout<<"ERROR in update_Hamiltonian_q_ethd: The entanglement option = "<<prms.entanglement_opt<<" is not avaialable\n";
    exit(0);
  }
  //cout<<"Stop here: \n"; exit(0);
}

void update_Hamiltonian_q_ethd(bp::dict prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params, MATRIX& invM){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_q_ethd(_prms, q, p, ham, py_funct, model_params, invM);

}

void update_Hamiltonian_q_ethd(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params){

  if(prms.entanglement_opt==0){    /* Nothing to do */   }
  else if(prms.entanglement_opt==1){   ham.add_ethd_adi(*dyn_var.q, *dyn_var.iM, 1);  }
  else if(prms.entanglement_opt==2){   ham.add_ethd3_adi(*dyn_var.q, *dyn_var.iM, prms.ETHD3_alpha, 1);  }
  else if(prms.entanglement_opt==22){  ham.add_ethd3_adi(*dyn_var.q, *dyn_var.p, *dyn_var.iM, prms.ETHD3_alpha, prms.ETHD3_beta, 1);  }
  else{
    cout<<"ERROR in update_Hamiltonian_q_ethd: The entanglement option = "<<prms.entanglement_opt<<" is not avaialable\n";
    exit(0);
  }

}

void update_Hamiltonian_q_ethd(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_q_ethd(_prms, dyn_var, ham, py_funct, model_params);
}



void update_nacs(dyn_control_params& prms, nHamiltonian& ham){
/**
  This function updates the internal (time-derivative, scalar) NACs in the hierarchy of 
  Hamiltonians
*/

  int isNBRA = prms.isNBRA;
  double dt = prms.dt;
  int nst = ham.nadi;
  int ntraj = ham.children.size();
  CMATRIX st(nst,nst);
  MATRIX st_re(nst, nst);
  MATRIX st_im(nst, nst);

  CMATRIX nac(nst, nst);
  MATRIX nac_re(nst, nst);
  MATRIX nac_im(nst, nst);



  if(isNBRA==1){
    if(prms.nac_update_method==2){

      st = ham.children[0]->get_time_overlap_adi();

      if(prms.nac_algo==0){        nac = 0.5*dt*(st-st.H());    }
      else if(prms.nac_algo==1){   
        st_re = st.real();
        nac_re = nac_npi(st_re, dt); 
        nac = CMATRIX(nac_re, nac_im);
      } 
 
      ham.children[0]->set_nac_adi_by_val(nac);
    }// if method == 2
        
  }// isNBRA == 1
  else{

    if(prms.nac_update_method==2){

      for(int traj=0; traj<ntraj; traj++){
        st = ham.children[traj]->get_time_overlap_adi();

        if(prms.nac_algo==0){        nac = 0.5*dt*(st-st.H());    }
        else if(prms.nac_algo==1){   
          st_re = st.real();
          nac_re = nac_npi(st_re, dt); 
          nac = CMATRIX(nac_re, nac_im);
        } 
 
        ham.children[traj]->set_nac_adi_by_val(nac);

      }// for traj
    }// for nac_update_method == 2

  }// else - isNBRA == 0

}



void update_Hamiltonian_p(dyn_control_params& prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM){

  /**
    Update of the vibronic Hamiltonian in response to changed p
  */

  // For the purpose of updating the NACs and Hvibs for just the quantum DOFs,
  // we'll reset the momenta for all other DOFs to zero, to effectively turn of
  // the effect of classical momenta on the NAC calculations (in case those derivative
  // couplings have been computed)
  int ndof, ntraj;
  vector<int>& which_dofs = prms.quantum_dofs;
  int n_active_dof = which_dofs.size();
  ndof = p.n_rows;
  ntraj = p.n_cols;

  MATRIX p_quantum_dof(ndof, ntraj);

  for(int idof = 0; idof < n_active_dof; idof++){
    int dof = which_dofs[idof];

    for(int itraj = 0; itraj < ntraj; itraj++){
      p_quantum_dof.set(dof, itraj,  p.get(dof, itraj) );
    }
  }

  //update_nacs(prms, ham);


  // Update NACs and Hvib for all trajectories
  if(prms.rep_tdse==0){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_dia(p_quantum_dof, invM, 0, 1);
    }
    ham.compute_hvib_dia(1);

  }
  else if(prms.rep_tdse==1){  

    if(prms.nac_update_method==0){ ;;  }
    else if(prms.nac_update_method==1){
      ham.compute_nac_adi(p_quantum_dof, invM, 0, 1); 
    }
    ham.compute_hvib_adi(1);
  }
}

void update_Hamiltonian_p(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham){

  update_Hamiltonian_p(prms, ham, *dyn_var.p, *dyn_var.iM);

}


void update_Hamiltonian_p(bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_p(_prms, ham, p, invM);

}


void update_Hamiltonian_p(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_p(_prms, ham, *dyn_var.p, *dyn_var.iM);

}




}// namespace libdyn
}// liblibra


