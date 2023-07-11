/*********************************************************************************
* Copyright (C) 2017-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Interactions.cpp
  \brief The file implements functions for molecular-mechanical Hamiltonian calculations as well
  as the classes for organizing such computations in an object-oriented way: basic functionality of the class
*/

#include "Interactions.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{


map<std::string, double> dict2map(boost::python::dict d){

  std::string key;
  double value;
  map<std::string, double> res;

  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    value = extract<double>(d.values()[i]);

    res.insert( std::pair<std::string, double>(key, value) ); 
  }

  return res;
}


//========================================================
// Basic functionality of the Interaction class
//========================================================

Interaction_N_Body::Interaction_N_Body(){  
/**
  Do nothing
*/
    Nbody = 0;

    is_active = 1;    // set it up to active by default
    int_type = -1;    // Undefined
    functional = -1;  // Undefined 

    energy = 0.0;     // Energy 
    stress_at = 0.0;  // stress tensor

    Hess = NULL;      // Hessian

}

Interaction_N_Body::Interaction_N_Body(int Nbody_){  
/**
  Allocates memory for the pointers
*/

    Nbody = Nbody_;

    r = vector<VECTOR*>(Nbody, NULL);
    t = vector<VECTOR*>(Nbody, NULL);
    f = vector<VECTOR*>(Nbody, NULL);
    q = vector<double*>(Nbody, NULL);


    is_active = 1;    // set it up to active by default
    int_type = -1;    // Undefined
    functional = -1;  // Undefined 

    energy = 0.0;     // Energy 
    stress_at = 0.0;  // stress tensor

    Hess = NULL;      // Hessian

}


void Interaction_N_Body::set_coords(VECTOR& r_, int indx){

    r[indx] = &r_; 

}// set_coords


void Interaction_N_Body::set_coords(vector<VECTOR*>& r_, vector<int>& indxs){

    address_subset<VECTOR>(r, r_, indxs);

}// set_coords

void Interaction_N_Body::set_coords(vector<VECTOR>& r_, vector<int>& indxs){

    address_subset<VECTOR>(r, r_, indxs);

}// set_coords



void Interaction_N_Body::set_transl(VECTOR& t_, int indx){

    t[indx] = &t_; 

}// set_transl

void Interaction_N_Body::set_transl(vector<VECTOR*>& t_, vector<int>& indxs){

    address_subset<VECTOR>(t, t_, indxs);

}// set_transl

void Interaction_N_Body::set_transl(vector<VECTOR>& t_, vector<int>& indxs){

    address_subset<VECTOR>(t, t_, indxs);

}// set_transl



void Interaction_N_Body::set_forces(VECTOR& f_, int indx){

    f[indx] = &f_; 

}// set_forces

void Interaction_N_Body::set_forces(vector<VECTOR*>& f_, vector<int>& indxs){

    address_subset<VECTOR>(f, f_, indxs);

}// set_forces

void Interaction_N_Body::set_forces(vector<VECTOR>& f_, vector<int>& indxs){

    address_subset<VECTOR>(f, f_, indxs);

}// set_forces


void Interaction_N_Body::set_charges(double& q_, int indx){

    q[indx] = &q_; 

}// set_charges


void Interaction_N_Body::set_charges(vector<double*>& q_, vector<int>& indxs){

    address_subset<double>(q, q_, indxs);

}// set_charges

void Interaction_N_Body::set_charges(vector<double>& q_, vector<int>& indxs){

    address_subset<double>(q, q_, indxs);

}// set_charges

void Interaction_N_Body::set_hessian(MATRIX& Hess_, vector<int>& Hess_stenc_){

    Hess = &Hess_;

    Hess_stenc = Hess_stenc_;

}// set Hessian


void Interaction_N_Body::set_functional(std::string f){
 ;;
}






}// namespace libhamiltonian_mm
}// namespace libatomistic 
}// liblibra
