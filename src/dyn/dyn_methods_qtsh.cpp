/*********************************************************************************
* Copyright (C) 2023 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_methods_qtsh.cpp
  \brief The file implements the quantum trajectory surface hopping (QTSH) method of Craig Martens

*/

#include "dyn_methods.h"
/*
#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "Dynamics.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"
#include "dyn_ham.h"

#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"
*/

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{


MATRIX compute_dkinemat(dyn_variables& dyn_var, nHamiltonian& ham){
/**
  This function computes the difference between the kinematic and canonical momenta

  p_{kin,j} = m dq_j/dt  = p_j - 2 \hbar * \sum_{n} \sum_{k<n}  d_{kn, j} \beta_{kn, j}

  So this function computes the p_{kin, j} - p_j

*/

  int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;
  int nadi = dyn_var.nadi;

  MATRIX dp(ndof, ntraj);

  for(int idof=0; idof<ndof; idof++){
    for(int itraj=0; itraj<ntraj; itraj++){

      double sum = 0.0;
      for(int n=0; n<nadi; n++){
        for(int k=0; k<n; k++){
          sum += ham.children[itraj]->dc1_adi[idof]->get(k, n).real() * dyn_var.dm_adi[itraj]->get(k, n).imag();
        }// for k
      }// for n

      dp.set(idof, itraj, -2.0*sum);

    }// for itraj
  }// for idof

  return dp;
}


}// libdyn
}// liblibra
