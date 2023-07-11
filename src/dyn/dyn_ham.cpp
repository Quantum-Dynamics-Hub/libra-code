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


void update_Hamiltonian_variables(dyn_control_params& prms, dyn_variables& dyn_var, 
                                  nHamiltonian& ham, nHamiltonian& ham_prev,
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

  int nadi = ham.nadi;
  int ntraj = ham.children.size();

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

    //============================== Time-overlaps ======================
    /// Don't update time-overlaps
    if(prms.time_overlap_method==0){  ;; }  // maybe it is already pre-computed and stored

    /// Compute the time-overlap directly, using previous MO vectors
    else if(prms.time_overlap_method==1){   // Explicitly compute it:

      CMATRIX st(nadi, nadi);
      for(int traj=0; traj<ntraj; traj++){
        st = ham_prev.children[traj]->get_basis_transform().H() * ham.children[traj]->get_basis_transform();
        ham.children[traj]->set_time_overlap_adi_by_val(st);
      }
    }// 1 

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
    else if(prms.nac_update_method==2 && update_type==0){

      // Don't update NACs and time-overlaps in response to p
      // if they are computed via the time-overlaps - only in response to 
      // change of q


      int isNBRA = prms.isNBRA;
      double dt = prms.dt;
      //int nst = ham.nadi;
      //int ntraj = ham.children.size();
      CMATRIX st(nadi,nadi);
      MATRIX st_re(nadi, nadi);
      MATRIX st_im(nadi, nadi);
 
      CMATRIX nac(nadi, nadi);
      MATRIX nac_re(nadi, nadi);
      MATRIX nac_im(nadi, nadi);

      for(int traj=0; traj<ntraj; traj++){
        st = ham.children[traj]->get_time_overlap_adi();

        if(prms.nac_algo==0){        nac = (0.5/dt)*(st-st.H());    }
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


  // Now apply the projectors to reflect any state switches that may have occurred in the meantime

  //cout<<"In update_Hamiltonian_variables... St = \n";
  //ham.children[0]->get_time_overlap_adi().show_matrix();

}


void update_Hamiltonian_variables(bp::dict prms, dyn_variables& dyn_var, 
                                  nHamiltonian& ham, nHamiltonian& ham_prev,
                                  bp::object py_funct, bp::object model_params, int update_type){

  dyn_control_params _prms;
  _prms.set_parameters(prms);

  update_Hamiltonian_variables(_prms, dyn_var, ham, ham_prev, py_funct, model_params, update_type);


}




}// namespace libdyn
}// liblibra


