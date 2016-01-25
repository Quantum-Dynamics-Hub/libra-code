/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

/**
  \file Dynamics_Ensemble.cpp
  \brief The file implements the functions for dynamics of ensemble of trajectories
    
*/

#include "Dynamics_Ensemble.h"

/// libdyn namespace
namespace libdyn{


void propagate_ensemble(double dt,Ensemble* ens,int opt){
/**
  \brief Propagate nuclear and electronic DOF for an ensemble of trajectories
  \param[in] dt Integration time step (evolution duration)
  \param[in,out] ens The pointer to the Ensemble object state of which is being propagated
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

*/

  // Ensemble of independent trajectories
  for(int i=0;i<ens->ntraj;i++){

    if(ens->is_active[i]){

      propagate_electronic(0.5*dt,&ens->el[i], ens->ham[i]);
      propagate_nuclear(dt, &ens->mol[i], &ens->el[i], ens->ham[i],opt);
      propagate_electronic(0.5*dt,&ens->el[i], ens->ham[i]);

    }
       
  }// for i

}

void propagate_ensemble(double dt,Ensemble& ens,int opt){
/**
  \brief Propagate nuclear and electronic DOF for an ensemble of trajectories - Python-friendly
  \param[in] dt Integration time step (evolution duration)
  \param[in,out] ens The Ensemble object state of which is being propagated
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

*/


  propagate_ensemble(dt, &ens, opt);
}






}// namespace libdyn


