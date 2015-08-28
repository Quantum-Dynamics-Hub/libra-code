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

#include "Surface_Hopping.h"
#include "Surface_Hopping_method1.h"
#include "Energy_and_Forces.h"

namespace libdyn{



void compute_hopping_probabilities_esh(Ensemble& ens, MATRIX* g, double dt, int use_boltz_factor,double T){

// esh - entangled surface hopping

// Assume ham internal processings are done outside
// assume matrix g is allocated


  const double kb = 3.166811429e-6; // Hartree/K
  int i,j,k,n, traj;
  double sum,g_ij,argg;


  // Verify if memory is allocated, if not generate
  if(g->num_of_elems!=ens.nelec * ens.nelec){cout<<"Matrix g is not allocated\n"; exit(0);  }


  // Now calculate the hopping probabilities
  // denominator
  vector<double> a_ii(ens.nelec, 0.0);

  for(i=0;i<ens.nelec;i++){
    for(traj=0; traj<ens.ntraj; traj++){

      a_ii[i] += ens.el[traj].q[i] * ens.el[traj].q[i] + ens.el[traj].p[i] * ens.el[traj].p[i]; // c_i^* * c_i

    }// for traj
    a_ii[i] /= ((double)ens.ntraj);

  }// for i

  // Kinetic energy spread:
  // Average values
  vector<double> p_ave(ens.nnucl,0.0);
  for(n=0; n<ens.nnucl; n++){  
    for(traj=0; traj<ens.ntraj; traj++){  

      p_ave[n] += ens.mol[traj].p[n];
    
    }// for traj
    p_ave[n] /= ((double)ens.ntraj);
  }// for n

  // Sigmas
  vector<double> p_sigma2(ens.nnucl,0.0);
  for(n=0; n<ens.nnucl; n++){  
    for(traj=0; traj<ens.ntraj; traj++){  

      p_sigma2[n] += (ens.mol[traj].p[n] - p_ave[n]) * (ens.mol[traj].p[n] - p_ave[n]);
    
    }// for traj
    p_sigma2[n] =  p_sigma2[n]/ ((double)ens.ntraj);
    
  }// for n

//  total_sigma = sqrt(0.1 * 2.0*2000.0*T*kb); // mass = 2000.0


  // nominator
  for(i=0;i<ens.nelec;i++){    
    double sum = 0.0;

    for(j=0;j<ens.nelec;j++){      
      if(i!=j){ // according to formula the diagonal probability P(i->i) should be zero
        double imHaij = 0.0;
        for(traj=0; traj<ens.ntraj; traj++){
 
          // Use very general expression:
          // Note: the sign here is not very obvious! Keep in mind:
          // Re(i*z) = -Im(z)  but  Im(i*z) = Re(z)
          // My formula is: P(i->j) = (2*dt/(hbar*|c_i|^2)) * ( -Im(H_ij * c_i^* * c_j) )
          double cij_re = ens.el[traj].q[i] * ens.el[traj].q[j] + ens.el[traj].p[i] * ens.el[traj].p[j];
          double cij_im = ens.el[traj].q[i] * ens.el[traj].p[j] - ens.el[traj].p[i] * ens.el[traj].q[j];
          imHaij += cij_re * ens.ham[traj]->Hvib(i,j).imag() + cij_im * ens.ham[traj]->Hvib(i,j).real(); // Im(H_ij * c_i^* * c_j)     

          // Coherence term:
          double coh = 0.0;
          for(n=0; n<ens.nnucl; n++){  
            coh += (2.0*ens.mol[traj].p[n] * ens.mol[traj].f[n]/p_sigma2[n]);
          }
          coh /= ( (ens.nelec-1) ); // equally distribute coherence over all states
          coh *= (ens.el[traj].q[i] * ens.el[traj].q[i] + ens.el[traj].p[i] * ens.el[traj].p[i]); 

          imHaij -= coh;

        }// for traj

        imHaij /= ((double)ens.ntraj);


        if(a_ii[i]<1e-8){ g_ij = 0.0; }  // avoid division by zero
        else{

          g_ij += -2.0*dt*imHaij/a_ii[i];  // This is a general case
          if(g_ij<0.0){  g_ij = 0.0; }

        }// else

        g->set(i,j,g_ij);
        sum = sum + g_ij;

      }// i!=j

      else{ g->set(i,j,0.0); }

    }// for j

    g->set(i,i,1.0 - sum);

  }// for i

}// compute probabilities_esh


}// namespace libdyn

