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


void update_Hamiltonian_variables(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                                  bp::object py_funct, bp::object model_params, int update_type){
/**
  This is a new function to replace all other "update_Hamiltonian functions"

  update_type:
     - 0: in response to changed q
     - 1: in response to changed p

  The key parameters in the `prms` are:

  time_overlap_method:
     - 0: based on the wavefunctions (the Hamiltonian shall have the basis_transform variables updated) 
     - 1: based on external calculations (the Hamiltonian shall have the time_overlap_adi member updated)

  ham_update_method:
     - 0: don't update any Hamiltonians
     - 1: recompute only diabatic Hamiltonian [ default ]
     - 2: recompute only adiabatic Hamiltonian

  ham_transform_method:
     - 0: don't do any transforms
     - 1: diabatic->adiabatic according to internal diagonalization [ default ]
     - 2: diabatic->adiabatic according to internally stored basis transformation matrix
     - 3: adiabatic->diabatic according to internally stored basis transformation matrix
     - 4: adiabatic->diabatic according to local diabatization method


  nac_update_method:
     - 0: don't update them (e.g. for simplest NAC)
     - 1: update according to changed momentum and existing derivative couplings [ default ]
     - 2: update according to time-overlaps (only time-derivative NACs)

  hvib_update_method:
     - 0: don't recompute Hvib (e.g. we read it)
     - 1: update according to changed NACs and energies
     
*/ 

  MATRIX& q = *dyn_var.q;
  MATRIX& p = *dyn_var.p;
  MATRIX& iM = *dyn_var.iM;


  if(update_type==0){
    //============================= Energies =========================
    // How to compute electronic Hamiltonian - this may already read a lot of 
    // other variables, such as time-overlaps, NAC, and Hvib
    if(prms.ham_update_method==0){ ;; }     
    else if(prms.ham_update_method==1){  
      ham.compute_diabatic(py_funct, q, model_params, 1);  
    }
    else if(prms.ham_update_method==2){  
      ham.compute_adiabatic(py_funct, q, model_params, 1);  
    }


    // Do the additional transformation between any reps, if needed
    if(prms.ham_transform_method==0){ ;; }
    else if(prms.ham_transform_method==1){
      ham.compute_adiabatic(1, 1);  // do internal diagonalization
    }
    else if(prms.ham_transform_method==2){
      // TO DO: Add dia->adi transformation according to the basis_transform matrix
    }
    else if(prms.ham_transform_method==3){
      // TO DO: Add adi->dia transformation according to the basis_transform matrix
    }
    else if(prms.ham_transform_method==4){
      // TO DO: Add adi->dia transformation according to local diabatization method
    }

        
    // Also add entanglement options - this is done in the adiabatic rep, so we'd need to compute 
    // all the adiabatic energies first
    if(prms.entanglement_opt==0){    /* Nothing to do */   }
    else if(prms.entanglement_opt==1){   ham.add_ethd_adi(q, iM, 1);  }
    else if(prms.entanglement_opt==2){   ham.add_ethd3_adi(q, iM, prms.ETHD3_alpha, 1);  }
    else if(prms.entanglement_opt==22){  ham.add_ethd3_adi(q, p, iM, prms.ETHD3_alpha, prms.ETHD3_beta, 1);  }
    else{
      cout<<"ERROR in update_Hamiltonian_q_ethd: The entanglement option = "<<prms.entanglement_opt<<" is not avaialable\n";
      exit(0);
    }
  }// update_type == 0


  if(update_type==0 || update_type==1){

    //========================== Couplings ===============================
    // Don't update NACs - they may have been read in step 1
    if(prms.nac_update_method==0){ ;;  }    

    // Compute NACs explicitly
    else if(prms.nac_update_method==1){
    
      // For the purpose of updating the NACs and Hvibs for just the quantum DOFs,
      // we'll reset the momenta for all other DOFs to zero, to effectively turn of
      // the effect of classical momenta on the NAC calculations (in case those derivative
      // couplings have been computed)
      int n_active_dof = prms.quantum_dofs.size();
      int ndof = dyn_var.ndof;
      int ntraj = dyn_var.ntraj;
      
      MATRIX p_quantum_dof(ndof, ntraj);
      for(auto dof: prms.quantum_dofs){
        for(int itraj = 0; itraj < ntraj; itraj++){  p_quantum_dof.set(dof, itraj,  p.get(dof, itraj) );  }
      }
    
      ham.compute_nac_dia(p_quantum_dof, iM, 0, 1);
      ham.compute_nac_adi(p_quantum_dof, iM, 0, 1);     
    }

    // Compute NAC from the time-overlaps
    else if(prms.nac_update_method==2){

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

    
    //========================== Vibronic Hamiltonian ===============================    
    // Don't update Hvib - perhaps because they are read from files in step 1
    if(prms.hvib_update_method==0){ ;;  }
    
    // Explicitly update Hvib
    else if(prms.hvib_update_method==1){    
      ham.compute_hvib_dia(1);
      ham.compute_hvib_adi(1);    
    }

  }// update_type==0 || update_type==1


}


void update_Hamiltonian_variables(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                                  bp::object py_funct, bp::object model_params, int update_type){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_variables(_prms, dyn_var, ham, py_funct, model_params, update_type);


}

void update_Hamiltonian_q(dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){
  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  int need_adiabatic_call = 0;

  // Propagation is in adiabatic basis but the Hamiltonian is in diabatic
  if(prms.rep_tdse==1 && prms.rep_ham==0){ need_adiabatic_call = 1; }

  // Using state-specific or Ehrenfest adiabatic forces, but the Hamiltonian is in diabatic rep
  if( (prms.force_method==1 || prms.force_method==2) && prms.rep_force == 1 && prms.rep_ham==0){ need_adiabatic_call = 1; }


  //------ Update the internals of the Hamiltonian object --------
  // We call the external function that would do the calculations

  if(prms.rep_ham==0){
    ham.compute_diabatic(py_funct, q, model_params, 1);
    ham.compute_adiabatic(1, 1);    
    //if(need_adiabatic_call){  ham.compute_adiabatic(1, 1);  }

  }
  else if(prms.rep_ham==1){
    ham.compute_adiabatic(py_funct, q, model_params, 1);
  }

//  ham.compute_hvib_dia(1);
//  ham.compute_hvib_adi(1);



/*
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
*/

}


void update_Hamiltonian_q(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params){

  /**
    Update of the vibronic Hamiltonian in response to changed q
  */

  update_Hamiltonian_q(prms, *dyn_var.q, ham, py_funct, model_params);

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


