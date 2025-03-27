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
#include "dyn_projectors.h"

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
     - 2: in response to changed electronic amplitudes

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
  int ndof = dyn_var.ndof;
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
//      cout<<" "<<q.n_cols<<"  "<<q.n_rows<<endl;
//      cout<<" "<<iM.n_cols<<"  "<<iM.n_rows<<endl;
//      cout<<" "<<p.n_cols<<"  "<<p.n_rows<<endl;
      //exit(0);
      ham.compute_diabatic(py_funct, q, model_params, 1);  
//      exit(0);
    }
    else if(prms.ham_update_method==2){  
      ham.compute_adiabatic(py_funct, q, model_params, 1);  
    }
//    exit(0);
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

    // Copy the transformation matrices to the dynamics variable
    dyn_var.update_basis_transform(ham); 


  }// update_type == 0

  //exit(0);


  if(update_type==0 || update_type==1){


    // Derivative NAC correction option:
    if(prms.do_nac_phase_correction==2){  // Experimental option to fix the phase of NACVs

      vector<int> traj_id(1,0);
      for(int traj=0; traj<ntraj; traj++){
        traj_id[0] = traj;
        CMATRIX Eadi(ham.children[traj]->get_ham_adi());
        MATRIX e_curr(ham.children[traj]->get_ham_adi().real());
        MATRIX e_prev(ham_prev.children[traj]->get_ham_adi().real());
        CMATRIX f_curr = ham.children[traj]->all_forces_adi(traj_id);
        CMATRIX f_prev = ham_prev.children[traj]->all_forces_adi(traj_id);
        //ham->get_hvib_adi().show_matrix();
        MATRIX pp(dyn_var.ndof, 1); 
        double dt = prms.dt;
        pp = p.col(traj);
        int act_state = dyn_var.act_states[traj];

        //vector<MATRIX> T_new(compute_F_cost_matrix_dof_resolved(f_curr, f_prev, e_curr, e_prev, pp, iM, dt, act_state));
        MATRIX T_new(compute_F_cost_matrix(f_curr, f_prev, e_curr, e_prev, pp, iM, dt, act_state).real() );

        for(int k=0; k<ndof;k++){
          for(int i=0; i<nadi;i++){
            for(int j=i+1; j<nadi; j++){
              //double sgn = T_new[k].get(i,j);
              double sgn = T_new.get(i,j);
              if(sgn > 0.5){ sgn = -1; } // change sign
              else{ sgn = 1.0; }
              
              ham.children[traj]->dc1_adi[k]->scale(i,j, sgn);
              ham.children[traj]->dc1_adi[k]->scale(j,i, sgn);
            }// for j
          }// for i

        }// for k
      }// for traj
    }// if correction

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

    if(prms.do_nac_phase_correction==1){  // Experimental option to fix the phase of NACs:
      for(int traj=0; traj<ntraj; traj++){
        for(int i=0; i<nadi;i++){
          for(int j=i+1; j<nadi; j++){
            double x1 = ham.children[traj]->nac_adi->get(i,j).real();
            double x2 = ham_prev.children[traj]->nac_adi->get(i,j).real();
            double sng =  SIGN(x1) * SIGN(x2); 
            ham.children[traj]->nac_adi->scale(i,j, sng);
            ham.children[traj]->nac_adi->scale(j,i, sng);
          }// for j
        }// for i
      }// for traj
    }// if correction

    
    //========================== Vibronic Hamiltonian ===============================    
    // Don't update Hvib - perhaps because they are read from files in step 1
    if(prms.hvib_update_method==0){ ;;  }
    
    // Explicitly update Hvib
    else if(prms.hvib_update_method==1){    
      ham.compute_hvib_dia(1);
      ham.compute_hvib_adi(1);    
    }

  }// update_type==0 || update_type==1


  if(update_type==2){  // in responce to updated electronic variables


  }// update_type==2


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


void update_time_overlaps(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham,  nHamiltonian& ham_prev){

  int nadi = ham.nadi;
  int ntraj = ham.children.size();

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

}

void update_proj_adi(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& Ham,  nHamiltonian& Ham_prev){
/**
  Just re-compute the proj_adi matrices
*/

  //======= Parameters of the dyn variables ==========
  int ntraj = dyn_var.ntraj;
  vector<int> traj_id(1,0);

  CMATRIX f_prev(dyn_var.nadi, dyn_var.ndof);
  CMATRIX f_curr(dyn_var.nadi, dyn_var.ndof);
  MATRIX iM(dyn_var.get_imass());
  MATRIX momenta(dyn_var.get_momenta());
  MATRIX p(dyn_var.ndof, 1);
  double dt = prms.dt;  

  double diff = 0.0;

  for(int itraj=0; itraj<ntraj; itraj++){
    int traj1 = itraj; // if(method >=100 && method <200){ traj1 = 0; }
    traj_id[0] = traj1;

    nHamiltonian* ham = Ham.children[traj1];
    nHamiltonian* ham_prev = Ham_prev.children[traj1];
    p = momenta.col(traj1);
    int act_state = dyn_var.act_states[traj1];

    //================= Basis re-expansion ===================
    CMATRIX P(ham->nadi, ham->nadi);
    CMATRIX T_new(*dyn_var.proj_adi[itraj]);

    P = ham->get_time_overlap_adi();  // (U_old.H() * U).H();

    if(prms.state_tracking_algo==-1){ // This is LD
      // More consistent with the new derivations:
      FullPivLU_inverse(P, T_new);   

      if( fabs( (P * T_new).tr().real() - P.n_cols) > 0.1 ){
        cout<<"Problem inverting time-overlap matrix\n";
        P.show_matrix("p_matrix.txt");
        exit(0);
      }
      T_new = orthogonalized_T( T_new );
    }
    else if(prms.state_tracking_algo==1 || prms.state_tracking_algo==2  || prms.state_tracking_algo==21 || 
            prms.state_tracking_algo==3 || prms.state_tracking_algo==32 || prms.state_tracking_algo==33){ // This is based on reordering + phase correction
      CMATRIX Eadi(ham->get_ham_adi());
      T_new = P;
      T_new = compute_projector(prms, Eadi, T_new);
    }
    else if(prms.state_tracking_algo==4){ // new experimental approach - based on forces
      //exit(0);
      CMATRIX Eadi(ham->get_ham_adi());
      MATRIX e_curr(ham->get_ham_adi().real());
      MATRIX e_prev(ham_prev->get_ham_adi().real());
      f_curr = ham->all_forces_adi(traj_id);
      f_prev = ham_prev->all_forces_adi(traj_id);
      //ham->get_hvib_adi().show_matrix();
      T_new = compute_F_cost_matrix(f_curr, f_prev, e_curr, e_prev, p, iM, dt, act_state);
      //T_new.show_matrix();
      //exit(0);
      T_new = compute_projector(prms, Eadi, T_new);  // CMATRIX compute_projector(dyn_control_params& prms, CMATRIX& Eadi, CMATRIX& St){
    //  T_new = orthogonalized_T( T_new );
    }

    *dyn_var.proj_adi[itraj] = T_new;

  }//for ntraj


}// reproject_basis

void update_proj_adi(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& Ham){
  update_proj_adi(prms, dyn_var, Ham, Ham);
}




}// namespace libdyn
}// liblibra


