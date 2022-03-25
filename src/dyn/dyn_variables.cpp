/*********************************************************************************
* Copyright (C) 2021 Alexey V. Akimov
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

dyn_variables::dyn_variables(int _ndia, int _nadi, int _ndof, int _ntraj){

/**

  This function initializes the default values of control parameters

*/

  cout<<"dyn_variables constructor!!!\n";  

  ///================= Dimension numbers =============
  ndia = _ndia;
  nadi = _nadi;
  ndof = _ndof;
  ntraj = _ntraj;

/*
  ///================= General variables, for OOP implementation ===================
  q = NULL;
  p = NULL;
  ampl_dia = NULL;
  ampl_adi = NULL;
  projectors = NULL;
*/

  ///================= A-FSSH ====================
  afssh_vars_status = 0;
  //dR = NULL;
  //dP = NULL;

  ///================= BCSH ====================
  bcsh_vars_status = 0;


}


void dyn_variables::allocate_afssh(){
//     cout<<"dyn_variables allocate_afssh!!!\n";  

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
//     cout<<"dyn_variables allocate_bcsh!!!\n";  

  if(bcsh_vars_status==0){

    reversal_events = new MATRIX(nadi, ntraj);
    bcsh_vars_status = 1;

  }

}// allocate_bcsh



dyn_variables::dyn_variables(const dyn_variables& x){     
     cout<<"dyn_variables copy constructor!!!\n";
     *this = x;
//    decoherence_rates = new MATRIX( *x.decoherence_rates );  

}



dyn_variables::~dyn_variables(){  
     cout<<"dyn_variables destructor!!!\n";
/*
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
*/

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

