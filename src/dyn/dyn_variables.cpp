/*********************************************************************************
* Copyright (C) 2021-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_variables.cpp
  \brief The file implements the methods to setup dynamical variable
*/

#include "dyn_variables.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;

void dyn_variables::allocate_electronic_vars(){

  if(electronic_vars_status==0){ 

    ampl_dia = new CMATRIX(ndia, ntraj);
    ampl_adi = new CMATRIX(nadi, ntraj);
    proj_adi = vector<CMATRIX*>(ntraj);
    dm_dia = vector<CMATRIX*>(ntraj);
    dm_adi = vector<CMATRIX*>(ntraj);

    for(int itraj=0; itraj<ntraj; itraj++){
      proj_adi[itraj] = new CMATRIX(nadi, nadi);
      proj_adi[itraj]->load_identity();
      dm_dia[itraj] = new CMATRIX(ndia, ndia);
      dm_adi[itraj] = new CMATRIX(nadi, nadi);
    }

    act_states = vector<int>(ntraj, 0);

    electronic_vars_status = 1;
  }

}

void dyn_variables::allocate_nuclear_vars(){

  if(nuclear_vars_status==0){ 

    iM = new MATRIX(ndof, 1);
    q = new MATRIX(ndof, ntraj);
    p = new MATRIX(ndof, ntraj);
    f = new MATRIX(ndof, ntraj);

    nuclear_vars_status = 1;
  }

}




dyn_variables::dyn_variables(int _ndia, int _nadi, int _ndof, int _ntraj){
/**

  This function initializes the default values of control parameters

*/

  ///================= Dimension numbers =============
  ndia = _ndia;
  nadi = _nadi;
  ndof = _ndof;
  ntraj = _ntraj;


  ///================= Electronic and nuclear variables, for OOP implementation ================
  electronic_vars_status = 0;
  allocate_electronic_vars();
    
  nuclear_vars_status = 0;
  allocate_nuclear_vars();

  ///================= A-FSSH ====================
  afssh_vars_status = 0;

  ///================= BCSH ====================
  bcsh_vars_status = 0;

  ///================= DISH ====================
  dish_vars_status = 0;

  ///================= FSSH2 ===================
  fssh2_vars_status = 0;
  
}


void dyn_variables::allocate_afssh(){


  if(afssh_vars_status==0){

    dR = vector< vector<CMATRIX*> >(ntraj, vector<CMATRIX*>(ndof, NULL) );
    dP = vector< vector<CMATRIX*> >(ntraj, vector<CMATRIX*>(ndof, NULL) );

    for(int itraj=0; itraj<ntraj; itraj++){
      for(int idof=0; idof<ndof; idof++){

        dR[itraj][idof] = new CMATRIX(nadi, nadi);
        dP[itraj][idof] = new CMATRIX(nadi, nadi);

      }
    }

    afssh_vars_status = 1;

  }

}// allocate_afssh



void dyn_variables::allocate_bcsh(){

  if(bcsh_vars_status==0){

    reversal_events = new MATRIX(nadi, ntraj);
    bcsh_vars_status = 1;

  }

}// allocate_bcsh


void dyn_variables::allocate_dish(){

  if(dish_vars_status==0){

    coherence_time = new MATRIX(nadi, ntraj);
    dish_vars_status = 1;

  }

}// allocate_dish


void dyn_variables::allocate_fssh2(){

  if(fssh2_vars_status==0){

    dm_dia_prev = vector<CMATRIX*>(ntraj);
    dm_adi_prev = vector<CMATRIX*>(ntraj);

    for(int itraj=0; itraj<ntraj; itraj++){
      dm_dia_prev[itraj] = new CMATRIX(ndia, ndia);
      dm_adi_prev[itraj] = new CMATRIX(nadi, nadi);
    }
    fssh2_vars_status = 1;
  }
}




dyn_variables::dyn_variables(const dyn_variables& x){     
  //cout<<"dyn_variables copy constructor\n";
  int itraj, idof;

  ndia = x.ndia;
  nadi = x.nadi;
  ndof = x.ndof;
  ntraj = x.ntraj;

  // copy content of electronic vars, only if initialized 
  if(x.electronic_vars_status==1){

    allocate_electronic_vars();

    *ampl_dia = *x.ampl_dia;
    *ampl_adi = *x.ampl_adi;
    for(itraj=0; itraj<ntraj; itraj++){
      *proj_adi[itraj] = *x.proj_adi[itraj];
      *dm_dia[itraj] = *x.dm_dia[itraj];
      *dm_adi[itraj] = *x.dm_adi[itraj];
    }
    act_states = x.act_states;

  }

  // copy content of nuclear vars, only if initialized 
  if(x.nuclear_vars_status==1){

    allocate_nuclear_vars();

    *iM = *x.iM;
    *q = *x.q;
    *p = *x.p;
    *f = *x.f;
  }

  // AFSSH vars - only if initialized
  if(x.afssh_vars_status==1){
    allocate_afssh();
    
    // Copy content
    for(itraj=0; itraj<ntraj; itraj++){
      for(idof=0; idof<ndof; idof++){
        *dR[itraj][idof] = *x.dR[itraj][idof];
        *dP[itraj][idof] = *x.dP[itraj][idof];
      }
    }

  }// if AFSSH vars

  // BCSH vars - only if initialized
  if(x.bcsh_vars_status==1){
    allocate_bcsh();

    // Copy content
    *reversal_events = *x.reversal_events;

  }// if BCSH vars

  // DISH vars - only if initialized
  if(x.dish_vars_status==1){
    allocate_dish();

    // Copy content
    *coherence_time = *x.coherence_time;

  }// if DISH vars

  // FSSH2 vars - only if initialized
  if(x.fssh2_vars_status==1){
    allocate_fssh2();
 
    // Copy content
    for(itraj=0; itraj<ntraj; itraj++){
      *dm_dia_prev[itraj] = *x.dm_dia_prev[itraj];
      *dm_adi_prev[itraj] = *x.dm_adi_prev[itraj];
    }    
  }// if FSSH2 vars
  
}// dyn_variables cctor



dyn_variables::~dyn_variables(){  
  //cout<<"dyn_variables destructor!!!\n";

  if(nuclear_vars_status==1){
    delete iM;
    delete q;
    delete p;
    delete f;
    nuclear_vars_status = 0;
  }

  if(electronic_vars_status==1){ 
    for(int itraj=0; itraj<ntraj; itraj++){
      delete proj_adi[itraj];
      delete dm_dia[itraj];
      delete dm_adi[itraj];
    }
    proj_adi.clear();
    dm_dia.clear();
    dm_adi.clear();
    
    delete ampl_dia;
    delete ampl_adi;

    act_states.clear();
    electronic_vars_status = 0;
  }

  if(afssh_vars_status==1){

    for(int itraj; itraj<ntraj; itraj++){
      for(int idof; idof<ndof; idof++){

        delete dR[itraj][idof];
        delete dP[itraj][idof];
      }// for idof

      dR[itraj].clear();
      dP[itraj].clear();

    }// for itraj

    dR.clear();
    dP.clear();

    afssh_vars_status = 0;

  }// AFSSH variables

  if(bcsh_vars_status==1){
    delete reversal_events;

    bcsh_vars_status = 0;
  }

  if(dish_vars_status==1){
    delete coherence_time;

    dish_vars_status = 0;
  }

  if(fssh2_vars_status==1){

    for(int itraj=0; itraj<ntraj; itraj++){
      delete dm_dia_prev[itraj];
      delete dm_adi_prev[itraj];
    }
    dm_dia_prev.clear();
    dm_adi_prev.clear();

    fssh2_vars_status = 0;
  }

}



CMATRIX dyn_variables::get_dm_adi(int i, int prev_steps){
  if(prev_steps==0){ return *dm_adi[i]; }
  else if(prev_steps==1){ return *dm_adi_prev[i]; }
  else{ ;; }
}

CMATRIX dyn_variables::get_dm_dia(int i, int prev_steps){
  if(prev_steps==0){ return *dm_dia[i]; }
  else if(prev_steps==1){ return *dm_dia_prev[i]; }
  else{ ;; }
}


void dyn_variables::set_parameters(bp::dict params){
/**
  Extract the parameters from the input dictionary
*/

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
//    if(key=="rep_tdse") { rep_tdse = bp::extract<int>(params.values()[i]); }
//    else if(key=="rep_ham") { rep_ham = bp::extract<int>(params.values()[i]);   }

  }
}


}// namespace libdyn
}// liblibra

