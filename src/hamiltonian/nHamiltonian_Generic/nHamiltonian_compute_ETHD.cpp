/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_ETHD.cpp
  \brief The file implements the calculations of the ETHD (Entangled Trajectories Hamiltonian Dynamics)
  terms 
    
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





double ETHD_energy(const MATRIX& q, const MATRIX& invM){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj;

  //============ Compute the averages =========
  MATRIX q_ave(ndof, 1);

  for(dof=0; dof<ndof; dof++){    
    for(traj=0; traj<ntraj; traj++){   q_ave.add(dof,0, q.get(dof, traj));  }
  }// for dof
  q_ave = q_ave / float(ntraj);

  //============ Compute the variances =========
  MATRIX s2(ndof, 1);

  for(dof=0; dof<ndof; dof++){    
    for(traj=0; traj<ntraj; traj++){   
      double s = q.get(dof, traj) - q_ave.get(dof,0);
      s2.add(dof,0, s*s);  
    }    
    s2.M[dof] = sqrt(s2.M[dof]/ float(ntraj));
  }// for dof


  //============ Compute the energy =========  
  double en = 0.0;

  for(dof=0; dof<ndof; dof++){    

    en += 0.125 * ntraj * invM.get(dof, 0) / s2.get(dof, 0);

  } 

  return en;


}


MATRIX ETHD_forces(const MATRIX& q, const MATRIX& invM){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  Returns:
  f - is a ndof x ntraj matrix that will contain the forces due to ETHD

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj;

  MATRIX f(ndof, ntraj);

  //============ Compute the averages =========
  MATRIX q_ave(ndof, 1);

  for(dof=0; dof<ndof; dof++){    
    for(traj=0; traj<ntraj; traj++){   q_ave.add(dof,0, q.get(dof, traj));  }
  }// for dof
  q_ave = q_ave / float(ntraj);

  //============ Compute the variances =========
  MATRIX s2(ndof, 1);

  for(dof=0; dof<ndof; dof++){    
    for(traj=0; traj<ntraj; traj++){   
      double s = q.get(dof, traj) - q_ave.get(dof,0);
      s2.add(dof,0, s*s);  
    }    
    s2.M[dof] = sqrt(s2.M[dof]/ float(ntraj));
  }// for dof


  //============ Compute the energy =========  
  for(dof=0; dof<ndof; dof++){    

    double s4 = s2.get(dof, 0); s4 = s4 * s4;
    double pref = 0.25 * invM.get(dof, 0) / s4;

    for(traj=0; traj<ntraj; traj++){
      
      f.set(dof, ntraj,  pref * (q.get(dof, traj) - q_ave.get(dof, traj) ) );

    }// for traj
  }//dof

  return f;

}







}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

