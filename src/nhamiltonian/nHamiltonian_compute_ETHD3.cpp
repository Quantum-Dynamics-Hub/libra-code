/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file nHamiltonian_compute_ETHD3.cpp
  \brief The file implements the calculations of the ETHD (Entangled Trajectories Hamiltonian Dynamics)
  terms. This new, experimental approach.
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif 

#include "nHamiltonian.h"
#include "../math_meigen/libmeigen.h"


/// liblibra namespace
namespace liblibra{

/// libnhamiltonian namespace 
namespace libnhamiltonian{


using namespace liblinalg;
using namespace libmeigen;





double ETHD3_energy(const MATRIX& q, const MATRIX& invM, double alp){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  alp - the coefficients of the Gaussian:  exp(-alp*(q_i-Q_i)^2)


  Complexity: O(Ntraj^2 x Ndof)

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj_k, traj_j;
  MATRIX dq(ndof, 1);

  //============ Compute the energy =========  
  double en = 0.0;

  for(traj_k=0; traj_k<ntraj; traj_k++){    


    double rho_k = 0.0; // total density at the position of trajectory k (Q_k)
    MATRIX d1rho_k(ndof, 1); // derivatives of the total density at the position of trajectory k (Q_k) w.r.t. q_{dof}
    MATRIX d2rho_k(ndof, 1); // second derivatives of the total density at the position of trajectory k (Q_k) w.r.t. q_{dof}
    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      dq = q.col(traj_k) - q.col(traj_j); 
      double g_kj = exp(- alp * (dq.T() * dq).get(0) );

      rho_k += g_kj;
      for(dof=0; dof<ndof; dof++){
        d1rho_k.add(dof, 0,  -2.0 * alp * dq.get(dof) * g_kj );
        d2rho_k.add(dof, 0,  (4.0 * alp * alp * dq.get(dof) * dq.get(dof) - 2.0 * alp) * g_kj );
      }// for dof
 
    }// for traj_j


    double A_k = 0.0;
    double B_k = 0.0;

    for(dof=0; dof<ndof; dof++){    

      A_k += d1rho_k.get(dof, 0) * invM.get(dof, 0) * d1rho_k.get(dof, 0);
      B_k += invM.get(dof, 0) * d2rho_k.get(dof, 0);
    }


    en += 0.125 * ( A_k / (rho_k * rho_k) -2.0 * B_k / rho_k); 

  }// for traj_k

  return en;

}



double ETHD3_energy(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  p - is a ndof x ntraj matrix of momenta
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  alp - the coefficients of the Gaussian:  exp(-alp*(q_i-Q_i)^2)
  bet - the coefficients of the Gaussian:  exp(-bet*(p_i-P_i)^2)


  Complexity: O(Ntraj^2 x Ndof)

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;

  int dof, traj_k, traj_j;
  MATRIX dq(ndof, 1);
  MATRIX dp(ndof, 1);

  //============ Compute the energy =========  
  double en = 0.0;

  for(traj_k=0; traj_k<ntraj; traj_k++){    


    double rho_k = 0.0; // total density at the position of trajectory k (Q_k)
    MATRIX d1rho_k(ndof, 1); // derivatives of the total density at the position of trajectory k (Q_k) w.r.t. q_{dof}
    MATRIX d2rho_k(ndof, 1); // second derivatives of the total density at the position of trajectory k (Q_k) w.r.t. q_{dof}
    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      dq = q.col(traj_k) - q.col(traj_j); 
      dp = p.col(traj_k) - p.col(traj_j); 
      double g_kj = exp(- alp * (dq.T() * dq).get(0) ) * exp(- bet * (dp.T() * dp).get(0) );

      rho_k += g_kj;
      for(dof=0; dof<ndof; dof++){
        d1rho_k.add(dof, 0,  -2.0 * alp * dq.get(dof) * g_kj );
        d2rho_k.add(dof, 0,  (4.0 * alp * alp * dq.get(dof) * dq.get(dof) - 2.0 * alp) * g_kj );
      }// for dof
 
    }// for traj_j


    double A_k = 0.0;
    double B_k = 0.0;

    for(dof=0; dof<ndof; dof++){    

      A_k += d1rho_k.get(dof, 0) * invM.get(dof, 0) * d1rho_k.get(dof, 0);
      B_k += invM.get(dof, 0) * d2rho_k.get(dof, 0);
    }


    en += 0.125 * ( A_k / (rho_k * rho_k) -2.0 * B_k / rho_k); 

  }// for traj_k

  return en;

}




MATRIX ETHD3_forces(const MATRIX& q, const MATRIX& invM, double alp){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  Returns:
  f - is a ndof x ntraj matrix that will contain the forces due to ETHD

  Complexity: O(Ntraj^2 x Ndof)

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int dof_i, dof_a, traj_k, traj_j, traj_n;

  MATRIX f(ndof, ntraj);
  MATRIX dq(ndof, 1);

  MATRIX g(ntraj, ntraj);
  MATRIX rho(1, ntraj);
  MATRIX A(1, ntraj);
  MATRIX B(1, ntraj);
  MATRIX d1rho(ndof, ntraj);
  MATRIX d2rho(ndof, ntraj);


  //============ Compute auxiliary variables =========  

  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      dq = q.col(traj_k) - q.col(traj_j); 
      double g_kj = exp(- alp * (dq.T() * dq).get(0) );

      g.set(traj_k, traj_j, g_kj);


      rho.add(0, traj_k, g_kj);

      for(dof_i=0; dof_i<ndof; dof_i++){
        d1rho.add(dof_i, traj_k,  -2.0 * alp * dq.get(dof_i) * g_kj );
        d2rho.add(dof_i, traj_k,  (4.0 * alp * alp * dq.get(dof_i) * dq.get(dof_i) - 2.0 * alp) * g_kj );
      }// for dof

    }// for traj_j
  }// for traj_k


  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      for(dof_i=0; dof_i<ndof; dof_i++){    

        A.add(0, traj_k, d1rho.get(dof_i, traj_j) * invM.get(dof_i, 0) * d1rho.get(dof_i, traj_j) );
        B.add(0, traj_k, invM.get(dof_i, 0) * d2rho.get(dof_i, traj_j));

      }// for dof_i
    }// for traj_j
  }// for traj_k




  // Complexity: O(Ntraj^3 x Ndof^2)

  for(traj_n=0; traj_n<ntraj; traj_n++){    
    for(dof_a=0; dof_a<ndof; dof_a++){    


      for(traj_k=0; traj_k<ntraj; traj_k++){    


        //========= Compute drho_k / dQ_an ===============
        double drhok_dQan = 0.0;

        for(traj_j=0; traj_j<ntraj; traj_j++){    
          if(traj_n==traj_j || traj_n==traj_k){ 

            double pref = (q.get(dof_a, traj_k) - q.get(dof_a, traj_j)) * g.get(traj_k, traj_j);
          
            if(traj_n==traj_k){   drhok_dQan += pref; }
            if(traj_n==traj_j){   drhok_dQan -= pref; }         

          }// if  
        }// traj_j

        drhok_dQan *= (-2.0*alp);

        // END: Compute drhok / dQ_an



        //========= Compute dA_k / dQ_an ===============
        double dAk_dQan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho'_ik / dQ_an ===============

          double drho_prime_ik_dQan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    
            if(traj_n==traj_j || traj_n==traj_k){ 

              double pref = -2.0*alp*(q.get(dof_i, traj_k) - q.get(dof_i, traj_j)) *  (q.get(dof_a, traj_k) - q.get(dof_a, traj_j)); 
              if(dof_i==dof_a){  pref += 1.0; }

              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){ drho_prime_ik_dQan += pref; }
              if(traj_n==traj_j){ drho_prime_ik_dQan -= pref; }

            }// if 
          }// traj_j

          drho_prime_ik_dQan *= (-2.0*alp);

          //==================================================
          
          dAk_dQan += 2.0*d1rho.get(dof_i, traj_k) * invM.get(dof_i, 0) * drho_prime_ik_dQan;

        }// dof_i
        // END: Compute dA_k / dQ_an


        //========= Compute dB_k / dQ_an ===============
        double dBk_dQan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho''_ik / dQ_an ===============

          double drho_dprime_ik_dQan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    
            if(traj_n==traj_j || traj_n==traj_k){ 

              double pref = 0.0;
              double dq_i_kj = (q.get(dof_i, traj_k) - q.get(dof_i, traj_j));
              double dq_a_kj = (q.get(dof_a, traj_k) - q.get(dof_a, traj_j));

              if(dof_i==dof_a){  pref = 2.0*dq_i_kj; }
              pref -= 2.0*alp* dq_a_kj  *dq_i_kj * dq_i_kj;
              pref += dq_a_kj;
              
              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){   drho_dprime_ik_dQan += pref; }
              if(traj_n==traj_j){   drho_dprime_ik_dQan -= pref; }

            }// if 
          }// traj_j

          drho_dprime_ik_dQan *= (4.0*alp*alp);

          //==================================================
          

          dBk_dQan += invM.get(dof_i, 0) * drho_dprime_ik_dQan;

        }// dof_i

        // END: Compute dB_k / dQ_an



        //======== Final formula =======================

        double Ak = A.get(0, traj_k);
        double Bk = B.get(0, traj_k);
        double rhok = rho.get(0, traj_k);
        double rhok2 = rhok * rhok;
        double rhok4 = rhok2 * rhok2;

        double val_k = (( dAk_dQan  + 2.0* Bk * drhok_dQan) * rhok2 - 2.0 * Ak * rhok * drhok_dQan  - 2.0 * dBk_dQan * rhok2 * rhok); 
        val_k /= rhok4;

        f.add(dof_a, traj_n, -0.125*val_k);

      }// for traj_k
    }// for dof_a
  }// for traj_n



  return f;

}



MATRIX ETHD3_forces(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet){
/**
  Compute the ETHD energy

  q - is a ndof x ntraj matrix of coordinates
  p - is a ndof x ntraj matrix of momenta
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  Returns:
  f - is a ndof x ntraj matrix that will contain the forces due to ETHD

  Complexity: O(Ntraj^2 x Ndof)

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int dof_i, dof_a, traj_k, traj_j, traj_n;

  MATRIX f(ndof, ntraj);
  MATRIX dq(ndof, 1);
  MATRIX dp(ndof, 1);

  MATRIX g(ntraj, ntraj);
  MATRIX rho(1, ntraj);
  MATRIX A(1, ntraj);
  MATRIX B(1, ntraj);
  MATRIX d1rho(ndof, ntraj);
  MATRIX d2rho(ndof, ntraj);


  //============ Compute auxiliary variables =========  

  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      dq = q.col(traj_k) - q.col(traj_j); 
      dp = p.col(traj_k) - p.col(traj_j); 
      double g_kj = exp(- alp * (dq.T() * dq).get(0) )  * exp(- bet * (dp.T() * dp).get(0) );

      g.set(traj_k, traj_j, g_kj);


      rho.add(0, traj_k, g_kj);

      for(dof_i=0; dof_i<ndof; dof_i++){
        d1rho.add(dof_i, traj_k,  -2.0 * alp * dq.get(dof_i) * g_kj );
        d2rho.add(dof_i, traj_k,  (4.0 * alp * alp * dq.get(dof_i) * dq.get(dof_i) - 2.0 * alp) * g_kj );
      }// for dof

    }// for traj_j
  }// for traj_k


  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      for(dof_i=0; dof_i<ndof; dof_i++){    

        A.add(0, traj_k, d1rho.get(dof_i, traj_j) * invM.get(dof_i, 0) * d1rho.get(dof_i, traj_j) );
        B.add(0, traj_k, invM.get(dof_i, 0) * d2rho.get(dof_i, traj_j));

      }// for dof_i
    }// for traj_j
  }// for traj_k




  // Complexity: O(Ntraj^3 x Ndof^2)

  for(traj_n=0; traj_n<ntraj; traj_n++){    
    for(dof_a=0; dof_a<ndof; dof_a++){    


      for(traj_k=0; traj_k<ntraj; traj_k++){    


        //========= Compute drho_k / dQ_an ===============
        double drhok_dQan = 0.0;

        for(traj_j=0; traj_j<ntraj; traj_j++){    
          if(traj_n==traj_j || traj_n==traj_k){ 

            double pref = (q.get(dof_a, traj_k) - q.get(dof_a, traj_j)) * g.get(traj_k, traj_j);
          
            if(traj_n==traj_k){   drhok_dQan += pref; }
            if(traj_n==traj_j){   drhok_dQan -= pref; }         

          }// if  
        }// traj_j

        drhok_dQan *= (-2.0*alp);

        // END: Compute drhok / dQ_an



        //========= Compute dA_k / dQ_an ===============
        double dAk_dQan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho'_ik / dQ_an ===============

          double drho_prime_ik_dQan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    
            if(traj_n==traj_j || traj_n==traj_k){ 

              double pref = -2.0*alp*(q.get(dof_i, traj_k) - q.get(dof_i, traj_j)) *  (q.get(dof_a, traj_k) - q.get(dof_a, traj_j)); 
              if(dof_i==dof_a){  pref += 1.0; }

              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){ drho_prime_ik_dQan += pref; }
              if(traj_n==traj_j){ drho_prime_ik_dQan -= pref; }

            }// if 
          }// traj_j

          drho_prime_ik_dQan *= (-2.0*alp);

          //==================================================
          
          dAk_dQan += 2.0*d1rho.get(dof_i, traj_k) * invM.get(dof_i, 0) * drho_prime_ik_dQan;

        }// dof_i
        // END: Compute dA_k / dQ_an


        //========= Compute dB_k / dQ_an ===============
        double dBk_dQan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho''_ik / dQ_an ===============

          double drho_dprime_ik_dQan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    
            if(traj_n==traj_j || traj_n==traj_k){ 

              double pref = 0.0;
              double dq_i_kj = (q.get(dof_i, traj_k) - q.get(dof_i, traj_j));
              double dq_a_kj = (q.get(dof_a, traj_k) - q.get(dof_a, traj_j));

              if(dof_i==dof_a){  pref = 2.0*dq_i_kj; }
              pref -= 2.0*alp* dq_a_kj  *dq_i_kj * dq_i_kj;
              pref += dq_a_kj;
              
              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){   drho_dprime_ik_dQan += pref; }
              if(traj_n==traj_j){   drho_dprime_ik_dQan -= pref; }

            }// if 
          }// traj_j

          drho_dprime_ik_dQan *= (4.0*alp*alp);

          //==================================================
          

          dBk_dQan += invM.get(dof_i, 0) * drho_dprime_ik_dQan;

        }// dof_i

        // END: Compute dB_k / dQ_an



        //======== Final formula =======================

        double Ak = A.get(0, traj_k);
        double Bk = B.get(0, traj_k);
        double rhok = rho.get(0, traj_k);
        double rhok2 = rhok * rhok;
        double rhok4 = rhok2 * rhok2;

        double val_k = (( dAk_dQan  + 2.0* Bk * drhok_dQan) * rhok2 - 2.0 * Ak * rhok * drhok_dQan  - 2.0 * dBk_dQan * rhok2 * rhok); 
        val_k /= rhok4;

        f.add(dof_a, traj_n, -0.125*val_k);

      }// for traj_k
    }// for dof_a
  }// for traj_n



  return f;

}




MATRIX ETHD3_friction(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet){
/**
  Compute the dE_{ETHD3} / dP_an

  q - is a ndof x ntraj matrix of coordinates
  p - is a ndof x ntraj matrix of momenta
  
  invM - is a ndof x 1 matrix of inverse masses of all DOFs

  Returns:
  f - is a ndof x ntraj matrix that will contain the forces due to ETHD

  Complexity: O(Ntraj^2 x Ndof)

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int dof_i, dof_a, traj_k, traj_j, traj_n;

  MATRIX f(ndof, ntraj);
  MATRIX dq(ndof, 1);
  MATRIX dp(ndof, 1);

  MATRIX g(ntraj, ntraj);
  MATRIX rho(1, ntraj);
  MATRIX A(1, ntraj);
  MATRIX B(1, ntraj);
  MATRIX d1rho(ndof, ntraj);
  MATRIX d2rho(ndof, ntraj);


  //============ Compute auxiliary variables =========  

  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      dq = q.col(traj_k) - q.col(traj_j); 
      dp = p.col(traj_k) - p.col(traj_j); 
      double g_kj = exp(- alp * (dq.T() * dq).get(0) )  * exp(- bet * (dp.T() * dp).get(0) );

      g.set(traj_k, traj_j, g_kj);


      rho.add(0, traj_k, g_kj);

      for(dof_i=0; dof_i<ndof; dof_i++){
        d1rho.add(dof_i, traj_k,  -2.0 * alp * dq.get(dof_i) * g_kj );
        d2rho.add(dof_i, traj_k,  (4.0 * alp * alp * dq.get(dof_i) * dq.get(dof_i) - 2.0 * alp) * g_kj );
      }// for dof

    }// for traj_j
  }// for traj_k


  // Complexity: O(Ntraj^2 x Ndof)

  for(traj_k=0; traj_k<ntraj; traj_k++){    
    for(traj_j=0; traj_j<ntraj; traj_j++){    

      for(dof_i=0; dof_i<ndof; dof_i++){    

        A.add(0, traj_k, d1rho.get(dof_i, traj_j) * invM.get(dof_i, 0) * d1rho.get(dof_i, traj_j) );
        B.add(0, traj_k, invM.get(dof_i, 0) * d2rho.get(dof_i, traj_j));

      }// for dof_i
    }// for traj_j
  }// for traj_k




  // Complexity: O(Ntraj^3 x Ndof^2)

  for(traj_n=0; traj_n<ntraj; traj_n++){    
    for(dof_a=0; dof_a<ndof; dof_a++){    


      for(traj_k=0; traj_k<ntraj; traj_k++){    


        //========= Compute drho_k / dP_an ===============
        double drhok_dPan = 0.0;

        for(traj_j=0; traj_j<ntraj; traj_j++){    
          if(traj_n==traj_j || traj_n==traj_k){ 

            double pref = (p.get(dof_a, traj_k) - p.get(dof_a, traj_j)) * g.get(traj_k, traj_j);
          
            if(traj_n==traj_k){   drhok_dPan += pref; }
            if(traj_n==traj_j){   drhok_dPan -= pref; }         

          }// if  
        }// traj_j

        drhok_dPan *= (-2.0*bet);

        // END: Compute drhok / dP_an



        //========= Compute dA_k / dP_an ===============
        double dAk_dPan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho'_ik / dP_an ===============

          double drho_prime_ik_dPan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    

            if(traj_n==traj_j || traj_n==traj_k){ 

              double pref = (q.get(dof_i, traj_k) - q.get(dof_i, traj_j)) *  (p.get(dof_a, traj_k) - p.get(dof_a, traj_j)); 

              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){ drho_prime_ik_dPan += pref; }
              if(traj_n==traj_j){ drho_prime_ik_dPan -= pref; }

            }// if 

          }// traj_j

          drho_prime_ik_dPan *= (4.0*alp * bet);

          //==================================================
          
          dAk_dPan += 2.0*d1rho.get(dof_i, traj_k) * invM.get(dof_i, 0) * drho_prime_ik_dPan;

        }// dof_i
        // END: Compute dA_k / dP_an


        //========= Compute dB_k / dQ_an ===============
        double dBk_dPan = 0.0;

        for(dof_i=0; dof_i<ndof; dof_i++){

          //========= Compute drho''_ik / dP_an ===============

          double drho_dprime_ik_dPan = 0.0;          

          for(traj_j=0; traj_j<ntraj; traj_j++){    

            if(traj_n==traj_j || traj_n==traj_k){ 

              
              double dq_i_kj = (q.get(dof_i, traj_k) - q.get(dof_i, traj_j));
              double dp_a_kj = (p.get(dof_a, traj_k) - p.get(dof_a, traj_j));

              double pref = (2.0*alp* dq_i_kj * dq_i_kj - 1.0);
              pref *= dp_a_kj;
              
              pref *= g.get(traj_k, traj_j);

              if(traj_n==traj_k){   drho_dprime_ik_dPan += pref; }
              if(traj_n==traj_j){   drho_dprime_ik_dPan -= pref; }

            }// if 

          }// traj_j

          drho_dprime_ik_dPan *= (-4.0*alp*bet);

          //==================================================
          

          dBk_dPan += invM.get(dof_i, 0) * drho_dprime_ik_dPan;

        }// dof_i

        // END: Compute dB_k / dP_an



        //======== Final formula =======================

        double Ak = A.get(0, traj_k);
        double Bk = B.get(0, traj_k);
        double rhok = rho.get(0, traj_k);
        double rhok2 = rhok * rhok;
        double rhok4 = rhok2 * rhok2;

        double val_k = (( dAk_dPan  + 2.0* Bk * drhok_dPan) * rhok2 - 2.0 * Ak * rhok * drhok_dPan  - 2.0 * dBk_dPan * rhok2 * rhok); 
        val_k /= rhok4;

        f.add(dof_a, traj_n, 0.125*val_k);

      }// for traj_k
    }// for dof_a
  }// for traj_n



  return f;

}





void nHamiltonian::add_ethd3_dia(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl){
/**
  Add ETHD3 energies and (optionally) forces to all the children Hamiltonians in the diabatic representation
*/

  complex<double> minus_one(-1.0, 0.0);

  if(der_lvl>=0){
    double en = ETHD3_energy(q, invM, alp);

    CMATRIX ethd_en(ndia, ndia);
    ethd_en.identity();
    ethd_en *= en;
    *ham_dia = ethd_en; 

    if(der_lvl>=1){
      MATRIX ethd_frcs(q.n_rows, q.n_cols);
      ethd_frcs = ETHD3_forces(q, invM, alp);

      for(int traj=0; traj<children.size(); traj++){

        for(int dof=0; dof<nnucl; dof++){

          for(int st=0; st<ndia; st++){

            children[traj]->d1ham_dia[dof]->add(st,st, ethd_frcs.get(dof, traj)*minus_one);

          }// for st
        }// dof
      }// traj

    }// der_lvl >=1

  }// der_lvl >=0

}


void nHamiltonian::add_ethd3_adi(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl){
/**
  Add ETHD3 energies and (optionally) forces to all the children Hamiltonians in the adiabatic representation
*/

  complex<double> minus_one(-1.0, 0.0);

  if(der_lvl>=0){
    CMATRIX ethd_en(nadi, nadi);
    ethd_en.identity();
    ethd_en *= ETHD3_energy(q, invM, alp);

    *ham_adi = ethd_en; 


    if(der_lvl>=1){
      MATRIX ethd_frcs(q.n_rows, q.n_cols);
      ethd_frcs = ETHD3_forces(q, invM, alp);

      for(int traj=0; traj<children.size(); traj++){

        for(int dof=0; dof<nnucl; dof++){

          for(int st=0; st<nadi; st++){

            children[traj]->d1ham_adi[dof]->add(st,st, ethd_frcs.get(dof, traj)*minus_one);

          }// for st
        }// dof
      }// traj

    }// der_lvl >=1

  }// der_lvl >=0

}



void nHamiltonian::add_ethd3_adi(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet, int der_lvl){
/**
  Add ETHD3 energies and (optionally) forces to all the children Hamiltonians in the adiabatic representation
*/

  complex<double> minus_one(-1.0, 0.0);

  if(der_lvl>=0){
    CMATRIX ethd_en(nadi, nadi);
    ethd_en.identity();
    ethd_en *= ETHD3_energy(q, p, invM, alp, bet);

    *ham_adi = ethd_en; 


    if(der_lvl>=1){
      MATRIX ethd_frcs(q.n_rows, q.n_cols);
      ethd_frcs = ETHD3_forces(q, p, invM, alp, bet);

      for(int traj=0; traj<children.size(); traj++){

        for(int dof=0; dof<nnucl; dof++){

          for(int st=0; st<nadi; st++){

            children[traj]->d1ham_adi[dof]->add(st,st, ethd_frcs.get(dof, traj)*minus_one);

          }// for st
        }// dof
      }// traj

    }// der_lvl >=1

  }// der_lvl >=0

}






}// namespace libnhamiltonian
}// liblibra

