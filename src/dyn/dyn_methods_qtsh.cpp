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


void update_forces_qtsh(dyn_variables& dyn_var, nHamiltonian& ham, dyn_control_params& prms){
  /**
    Add the nonclassical contribution from the kinematic mometum in QTSH
  */

  int ntraj = dyn_var.ntraj;
  int ndof = dyn_var.ndof;
  int nst = dyn_var.nadi;
  MATRIX& invM = *dyn_var.iM;
  
  CMATRIX C(nst, 1);
  CMATRIX Coeff(nst, ntraj);
  
  MATRIX f_nc(ndof,ntraj);

  Coeff = *dyn_var.ampl_adi;
  
  for(int traj=0; traj<ntraj; traj++){
    C = Coeff.col(traj);

    CMATRIX ham_adi(nst, nst);
    ham_adi = ham.children[traj]->get_ham_adi();
    
    CMATRIX dc1_adi(nst, nst);
    for(int idof=0; idof<ndof; idof++){
      dc1_adi = ham.children[traj]->get_dc1_adi(idof);
      for(int i=0; i<nst; i++){
        for(int j=0; j<nst; j++){
          if(i<=j){continue;}
          f_nc.add(idof, traj, 2.0*dc1_adi.get(i,j).real() * (ham_adi.get(i,i) - ham_adi.get(j,j)).real() 
            * dyn_var.dm_adi[traj]->get(i,j).real() );
        }
      }
    }
  }//traj
  
  *dyn_var.f += f_nc;

}

void update_deco_factor_qtsh(dyn_variables& dyn_var, nHamiltonian& ham, dyn_control_params& prms){
/**
  This function computes the global decoherence factor used to adjust the hopping probability in QTSH
*/
  int ndof = dyn_var.ndof;
  int ntraj = dyn_var.ntraj;
  int nadi = dyn_var.nadi;
  
  double dt = prms.dt;

  // Compute the proxy phase of each trajectory
  for(int traj=0; traj<ntraj; traj++){
    CMATRIX ham_adi(nadi, nadi);
    ham_adi = ham.children[traj]->get_ham_adi();
    
    for(int i=0; i<nadi; i++){
      for(int j=0; j<nadi; j++){
        double omega = ham_adi.get(i,i).real() - ham_adi.get(j,j).real();
        dyn_var.qtsh_proxy_phase[traj]->add(i,j, omega*dt);
      }
    }
  }
  
  // Compute the global decoherence factor
  MATRIX avg_ph(nadi, nadi); 
  MATRIX avg_ph2(nadi, nadi);
  MATRIX temp(nadi, nadi);
  
  for(int traj=0; traj<ntraj; traj++){
    MATRIX phi(*dyn_var.qtsh_proxy_phase[traj]);
    avg_ph += phi;
    
    temp.dot_product(phi, phi);
    avg_ph2 += temp;
  }
  avg_ph /= ntraj;
  avg_ph2 /= ntraj;

  for(int i=0; i<nadi; i++){
    for(int j=0; j<nadi; j++){
      dyn_var.qtsh_deco_factor->set(i,j, exp( -0.5*(avg_ph2.get(i,j) - pow(avg_ph.get(i,j), 2.0)) ) );
    }
  }

}

}// libdyn
}// liblibra
