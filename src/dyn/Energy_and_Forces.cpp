/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Energy_and_Forces.cpp
  \brief The file implements the functions for dynamics-immediate energy calculations

  The "dynamics-immediate" means the energies and forces computed and organized to 
  be used in dynamical calculations of different type - classical, quantum, quantum-classical.
    
*/

#include "Energy_and_Forces.h"
#include "dyn_projectors.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace 
namespace libdyn{


double compute_kinetic_energy(MATRIX& p, MATRIX& invM, vector<int>& which_dofs){
/**
  \brief Compute kinetic energy of nuclear DOFs
  \param[in] p [ndof x ntraj] Momenta of ntraj replicas of the system with ndof 
  nuclear DOFs
  \param[in] invM [ndof x 1] Matrix of inverted masses of all DOFs

  In the case of more than 1 trajectories supplied into this function, the 
  total kinetic energy of all trajectories is returned.

  This is the classical nuclear kinetic energy
*/
  int ndof = which_dofs.size(); ///p.n_rows;
  int ntraj = p.n_cols;

  double Ekin = 0.0;
  
    
  for(int idof=0; idof < ndof; idof++){ 

    int dof = which_dofs[idof];

    double sum = 0.0;
    for(int traj=0; traj < ntraj; traj++){    
      sum = p.get(dof, traj) * p.get(dof, traj);
    }
    Ekin += sum * invM.get(dof, 0);
  }
  Ekin *= 0.5;

  return Ekin;

}

double compute_kinetic_energy(MATRIX& p, MATRIX& invM){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof);
  for(int i = 0; i < ndof; i++){  which_dofs[i] = i; }

  return compute_kinetic_energy(p, invM, which_dofs);

}


vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM, vector<int>& which_dofs){
/**
  \brief Compute trajectory-resolved kinetic energies of nuclear DOFs

  \param[in] p [ndof x ntraj] Momenta of ntraj replicas of the system with ndof 
  nuclear DOFs
  \param[in] invM [ndof x 1] Matrix of inverted masses of all DOFs

  This is the classical nuclear kinetic energy
*/
  int ndof = which_dofs.size(); ///p.n_rows;
  int ntraj = p.n_cols;

  vector<double> Ekin(ntraj, 0.0);
  
  for(int traj=0; traj < ntraj; traj++){    

    double sum = 0.0;

    for(int idof=0; idof < ndof; idof++){ 

      int dof = which_dofs[idof];
      sum += p.get(dof, traj) * invM.get(dof, 0) * p.get(dof, traj);
    }
    sum *= 0.5;

    Ekin[traj] = sum;

  }// for traj

  return Ekin;

}



vector<double> compute_kinetic_energies(MATRIX& p, MATRIX& invM){

  int ndof = p.n_rows;
  vector<int> which_dofs(ndof);
  for(int i = 0; i < ndof; i++){  which_dofs[i] = i; }

  return compute_kinetic_energies(p, invM, which_dofs);

}





CMATRIX tsh_indx2ampl(vector<int>& res, int nstates){
/**
  Convert the active state index to the vector with occupation numbers
  This is done for 1 or many trajectories
*/

  int ntraj = res.size();
  CMATRIX ampl(nstates, ntraj);

  
  for(int traj = 0; traj < ntraj; traj++){

    ampl.set(res[traj], traj, complex<double>(1.0, 0.0));

  }

  return ampl;

}


double average_potential_energy(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){

  vector<double> res(dyn_vars.ntraj, 0.0);
  res = potential_energies(prms, dyn_vars, ham);

  double ave = 0.0;
  for(int itraj = 0; itraj < dyn_vars.ntraj; itraj++){  ave += res[itraj]; }
  ave /= double(dyn_vars.ntraj);
  
  return ave;
}

double average_potential_energy(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  return average_potential_energy(prms, dyn_vars, ham);
}


vector<double> potential_energies(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){

  int itraj;
  int ntraj = dyn_vars.ntraj;
  vector<int> id(2,0);

  vector<double> res(ntraj, 0.0);

  if(prms.force_method==0){  // No forces

    // Don't compute forces at all - e.g. in NBRA
    for(itraj=0; itraj<ntraj; itraj++){ res[itraj] = 0.0; }

  }// NBRA

  else if(prms.force_method==1){   // State-specific forces

    // TSH or adiabatic (including excited states)
    // state-specific forces

    vector<int> effective_states(dyn_vars.act_states);
    //vector<int> eff_states(dyn_vars.act_states);

    if(prms.enforce_state_following==1){ 
      for(itraj=0; itraj<ntraj; itraj++){ effective_states[itraj] = prms.enforced_state_index; }
    }

    //eff_states = update_active_states(effective_states, dyn_vars.proj_adi);

    // Diabatic 
    if(prms.rep_force==0){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        int ist = effective_states[itraj];
        res[itraj] = ham.get_ham_dia(id).get(ist, ist).real();   
      }
    }
    // Adiabatic 
    else if(prms.rep_force==1){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        int ist = effective_states[itraj];

        //CMATRIX& T = *dyn_vars.proj_adi[itraj];
        int nst = dyn_vars.nadi;
/*
        CMATRIX T(nst, nst);
        T = orthogonalized_T( *dyn_vars.proj_adi[itraj] );
        res[itraj] = (T.H() * ham.get_ham_adi(id) * T).get(ist, ist).real();   
*/
        res[itraj] = ham.get_ham_adi(id).get(ist,ist).real();

//        if(dyn_vars.q->get(0,0)>-1.0 && dyn_vars.q->get(0,0)<3.0 ){
//          cout<<"E_pot = "<<ham.get_ham_adi(id).get(ist, ist).real()<<"  "<<res[itraj]<<endl;
//          cout<<" states: active = "<<effective_states[0]<<"  corrected eff states = "<<eff_states[0]<<endl;
//        }

      }
    }

  }// TSH && adiabatic

  else if(prms.force_method==2){  // Ehrenfest forces

    // Diabatic 
    if(prms.rep_force==0){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        res[itraj] = ham.Ehrenfest_energy_dia(*dyn_vars.ampl_dia, id).real();
      }
    }
    // Adiabatic 
    else if(prms.rep_force==1){  
      for(itraj=0; itraj<ntraj; itraj++){
        id[1] = itraj;
        res[itraj] = ham.Ehrenfest_energy_adi(*dyn_vars.ampl_adi, id).real();
      }
    }
  
  }// Ehrenfest


  return res;
}

vector<double> potential_energies(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  return potential_energies(prms, dyn_vars, ham);

}

//MATRIX aux_get_forces(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){
void update_forces(dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham){
  /**
    Compute the force depending on the method used
  */

  //cout<<"aux_get_forces \n";

  int ndof = ham.nnucl;
  int nst = ham.nadi;
  int ntraj = dyn_vars.ntraj;

  //MATRIX F(ndof, ntraj);
  CMATRIX f_all(nst, ndof);
  CMATRIX f_diag(nst, nst);
  CMATRIX f(ndof, 1);
  vector<int> id(2, 0); 
  

  if(prms.force_method==0){  // No forces

    // Don't compute forces at all - e.g. in NBRA

  }// NBRA

  else if(prms.force_method==1){   // State-specific forces

    // TSH or adiabatic (including excited states)
    // state-specific forces   
    vector<int> effective_states(dyn_vars.act_states);
//    vector<int> eff_states(dyn_vars.act_states);
  
    if(prms.enforce_state_following==1){ // NBRA-like enforcement: adiabatic dynamics, in terms of forces 
      for(int itraj=0; itraj<ntraj; itraj++){ effective_states[itraj] = prms.enforced_state_index; }
    }

    //eff_states = update_active_states(effective_states, dyn_vars.proj_adi);

    // Diabatic 
    if(prms.rep_force==0){  *dyn_vars.f = ham.forces_dia(effective_states).real();   }

    // Adiabatic 
    else if(prms.rep_force==1){  
      *dyn_vars.f = ham.forces_adi(effective_states).real(); //- older approach, without reordering 


//        if(dyn_vars.q->get(0,0)>-1.0 && dyn_vars.q->get(0,0)<3.0 ){
//        cout<<"q = \n"; dyn_vars.q->show_matrix();  dyn_vars.proj_adi[0]->show_matrix(); 
//        cout<<"f = \n"; dyn_vars.f->show_matrix();
//        cout<<" states: active = "<<effective_states[0]<<endl;
        //cout<<"==\n";
//        }


/*
      for(int traj=0; traj<ntraj; traj++){

        id[0] = 0; id[1] = traj;
        f_all = ham.all_forces_adi( id ); // forces on all adiabatic states of traj `traj`
        //CMATRIX& T = *dyn_vars.proj_adi[traj];
        
        CMATRIX T(nst, nst);       
        T = orthogonalized_T( *dyn_vars.proj_adi[traj] );

        if(dyn_vars.q->get(0,0)>-1.0 && dyn_vars.q->get(0,0)<3.0 ){
        cout<<"q = \n"; dyn_vars.q->show_matrix();  dyn_vars.proj_adi[traj]->show_matrix(); T.show_matrix();
        cout<<"initial f = \n"; dyn_vars.f->show_matrix();
        //cout<<"==\n";
        }

        for(int idof=0; idof<ndof; idof++){  
          for(int ist=0; ist<nst; ist++){  f_diag.set(ist, ist, f_all.get(ist, idof) ); }
          f_diag =  T.H() * f_diag * T;
          dyn_vars.f->set(idof, traj, f_diag.get( effective_states[traj],  effective_states[traj] ).real() );  
        }

        if(dyn_vars.q->get(0,0)>-1.0 && dyn_vars.q->get(0,0)<3.0 ){
          cout<<"final f = \n"; dyn_vars.f->show_matrix();
        }

      } // for traj
*/

    }// rep_force == 1
  }// TSH && adiabatic

  else if(prms.force_method==2){  // Ehrenfest forces
    // Diabatic 
    if(prms.rep_force==0){  *dyn_vars.f = ham.Ehrenfest_forces_dia(*dyn_vars.ampl_dia, 1).real();   }

    // Adiabatic 
    //cout<<"Ampl = \n";
    //dyn_vars.ampl_adi->show_matrix();
    else if(prms.rep_force==1){ *dyn_vars.f = ham.Ehrenfest_forces_adi(*dyn_vars.ampl_adi, dyn_vars.proj_adi, 1).real(); }

    //cout<<"Ampl = "<<dyn_vars.ampl_adi<<endl;
    
    //dyn_vars.ampl_adi->show_matrix();

    //cout<<"Ehrenfest force \n"; F.show_matrix();
  
  }// Ehrenfest


//  return F;
}

//MATRIX aux_get_forces(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){
void update_forces(bp::dict params, dyn_variables& dyn_vars, nHamiltonian& ham){

  dyn_control_params prms;
  prms.set_parameters(params);

  update_forces(prms, dyn_vars, ham);
//  return aux_get_forces(prms, dyn_vars, ham);
  
}



vector<CMATRIX> get_Eadi(nHamiltonian& ham){

  int nst = ham.nadi;
  int ntraj = ham.children.size();

  vector<CMATRIX> Eadi(ntraj, CMATRIX(nst, nst));

  for(int traj=0; traj<ntraj; traj++){
    Eadi[traj] = ham.children[traj]->get_ham_adi(); 
  }

  return Eadi;

}

vector<CMATRIX> get_Eadi(nHamiltonian* ham){

  int nst = ham->nadi;
  int ntraj = ham->children.size();

  vector<CMATRIX> Eadi(ntraj, CMATRIX(nst, nst));

  for(int traj=0; traj<ntraj; traj++){
    Eadi[traj] = ham->children[traj]->get_ham_adi();
  }

  return Eadi;

}



}// namespace libdyn
}// liblibra
