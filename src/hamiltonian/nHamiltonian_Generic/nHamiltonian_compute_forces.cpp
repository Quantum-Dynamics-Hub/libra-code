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
  \file nHamiltonian_compute_forces.cpp
  \brief The file implements computations of various types of forces and force tensors
    
*/


#include <stdlib.h>

#include "nHamiltonian.h"
#include "../../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{

/// libhamiltonian namespace 
namespace libhamiltonian{

/// libhamiltonian_generic namespace 
namespace libhamiltonian_generic{

using namespace liblinalg;
using namespace libmeigen;



CMATRIX nHamiltonian::forces_adi(vector<int>& act_states){
/**
  This only applies for a nHamiltonian object with a set of children whose size is 
  equal to the number of integers in `act_states` input

  Return: -d1ham_adi shaped as CMATRIX(ndof, ntraj) - for `act_states` states
*/
  int ntraj = act_states.size();   
  CMATRIX res(nnucl, ntraj);
 
  if(children.size()==ntraj){

    for(int traj=0; traj<ntraj; traj++){
      int st = act_states[traj];

      for(int dof=0; dof<nnucl; dof++){
        res.set(dof, traj, -children[traj]->d1ham_adi[dof]->get(st, st) );
      }// for dof 
    }// for traj
  }
  else{
    cout<<"ERROR: the size of the input is different from the number of children\n"; exit(0); 
  }

  return res;
}



CMATRIX nHamiltonian::forces_dia(vector<int>& act_states){
/**
  This only applies for a nHamiltonian object with a set of children whose size is 
  equal to the number of integers in `act_states` input

  Return: -d1ham_dia shaped as CMATRIX(ndof, ntraj) - for `act_states` states
*/

  int ntraj = act_states.size();   
  CMATRIX res(nnucl, ntraj);
 
  if(children.size()==ntraj){

    for(int traj=0; traj<ntraj; traj++){
      int st = act_states[traj];

      for(int dof=0; dof<nnucl; dof++){
        res.set(dof, traj, -children[traj]->d1ham_dia[dof]->get(st, st) );
      }// for dof 
    }// for traj
  }
  else{
    cout<<"ERROR: the size of the input is different from the number of children\n"; exit(0); 
  }

  return res;
}




}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

