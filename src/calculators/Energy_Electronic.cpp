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

#include "Energy_Electronic.h"

namespace libcalculators{
 
double energy_elec(MATRIX* Pao,MATRIX* Hao,MATRIX* Fao){
// Compute electronic energy (true for HF-derived methods: HF, CNDO, CNDO/2, INDO)
// this general formula is also true for EHT (F = Hcore, so the energy is simply a weighted sum of the eigenvalues)

  int Norb = Pao->num_of_cols;
  
  int i,ii,j;
  double Eelec = 0.0;

  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){

      Eelec += Pao->M[i*Norb+j]*(Hao->M[i*Norb+j] + Fao->M[i*Norb+j]);

    }// for i
  }// for j

  Eelec *= 0.5;

  return Eelec;

}// double energy_elec(...)

double energy_elec(MATRIX Pao,MATRIX Hao,MATRIX Fao){

  return energy_elec(&Pao,&Hao,&Fao);
}


double energy_elec(MATRIX* P_alp, MATRIX* P_bet, 
                   MATRIX* Hao_alp, MATRIX* Hao_bet,
                   MATRIX* Fao_alp, MATRIX* Fao_bet,
                   MATRIX* dFao_alp_dP_alp, MATRIX* dFao_alp_dP_bet,
                   MATRIX* dFao_bet_dP_alp, MATRIX* dFao_bet_dP_bet,
                   MATRIX* temp
                  ){
// Compute electronic energy (true for HF-derived methods: HF, CNDO, CNDO/2, INDO)
// this general formula is also true for EHT (F = Hcore, so the energy is simply a weighted sum of the eigenvalues)
// temp - is just a temporary array - preallocate it before calling this function
//        this will give some acceleration if the energy function is called very often
// Scaling - O(N^3)

  int Norb = P_alp->num_of_cols;
  
  int i,ii,j;
  double Eelec_alp = 0.0;
  double Eelec_bet = 0.0;
  double Eelec_tot = 0.0;

  //-------------------- Handle alpha part ---------------------------------------
  *temp = 0.0; 
  if(dFao_alp_dP_alp->max_elt()>1e-10){   *temp += *P_alp * *dFao_alp_dP_alp;    }
  if(dFao_alp_dP_bet->max_elt()>1e-10){   *temp += *P_bet * *dFao_alp_dP_bet;    }
   
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      Eelec_alp += P_alp->M[i*Norb+j]*(Hao_alp->M[i*Norb+j] +( Fao_alp->M[i*Norb+j] + temp->M[i*Norb+j]));
    }// for i
  }// for j
  Eelec_alp *= 0.5;


  //-------------------- Handle beta part ---------------------------------------
  *temp = 0.0; 
  if(dFao_bet_dP_alp->max_elt()>1e-10){   *temp += *P_alp * *dFao_bet_dP_alp;    }
  if(dFao_bet_dP_bet->max_elt()>1e-10){   *temp += *P_bet * *dFao_bet_dP_bet;    }
   
  for(i=0;i<Norb;i++){
    for(j=0;j<Norb;j++){
      Eelec_bet += P_bet->M[i*Norb+j]*(Hao_bet->M[i*Norb+j] +( Fao_bet->M[i*Norb+j] + temp->M[i*Norb+j]));
    }// for i
  }// for j
  Eelec_bet *= 0.5;


  //------------------------- Total ----------------------------------------------
  Eelec_tot = Eelec_alp + Eelec_bet;


  return Eelec_tot;

}// double energy_elec(...)

double energy_elec(MATRIX P_alp, MATRIX P_bet, 
                   MATRIX Hao_alp, MATRIX Hao_bet,
                   MATRIX Fao_alp, MATRIX Fao_bet,
                   MATRIX dFao_alp_dP_alp, MATRIX dFao_alp_dP_bet,
                   MATRIX dFao_bet_dP_alp, MATRIX dFao_bet_dP_bet,
                   MATRIX temp
                  ){

  return energy_elec(&P_alp,&P_bet,&Hao_alp,&Hao_bet,&Fao_alp,&Fao_bet,
                     &dFao_alp_dP_alp, &dFao_alp_dP_bet, &dFao_bet_dP_alp, &dFao_bet_dP_bet, &temp);
}


}// namespace libcalculators

