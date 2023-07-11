/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Potentials_vdw.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{


double Vdw_LJ(VECTOR& ri,VECTOR& rj,          /*Inputs*/
              VECTOR& fi,VECTOR& fj,          /*Outputs*/
              double sigma, double epsilon){  /*Parameters*/
//****************** double Lennard-Jones potential **************************
//*                                                                          *
//*       E = epsilon_min*[(sigma_min/r_ij)^12 - 2(sigma_min/r_ij)^6]        *
//*                                                                          *
//* E(sigma_min) = min E = -epsilon_min                                      *
//*                 r                                                        *
//* Another version of LJ potential is:                                      *
//*       E = epsilon_zero*[(sigma_zero/r_ij)^12 - (sigma_zero/r_ij)^6]      *
//* E(sigma_zero) = 0                                                        *
//* But since the latter form is related to the former one by change of      *
//* parameters it is enough to implement only the former one:                *
//* Parameters are related:                                                  *
//* epsilon_zero = 4*epsilon_min                                             *
//* sigma_zero   = sigma_min/pow(sigma_min, (1/6))                           *
//****************************************************************************
  double energy,r2,r6,r12,d2;
  VECTOR rij = ri - rj;
  d2 = rij.length2();
  r2 = (sigma*sigma/d2);
  r6 = r2*r2*r2;
  r12= r6*r6;
  energy = epsilon*(r12-2.0*r6);
  fi = 12.0*epsilon*(r12 - r6)*(rij/d2);
  fj = -fi;
  return energy;
}

double Vdw_Buffered14_7(VECTOR& ri,VECTOR& rj,          /*Inputs*/
                        VECTOR& fi,VECTOR& fj,          /*Outputs*/
                        double sigma, double epsilon){  /*Parameters*/
//****************** double Buffered 14-7 potential *************************
//*
//*  u = epsilon * ( 1.07*sigma_{ij}/(r_ij + 0.07 sigma_ij) )^7 [ (1.12 sigma_ij^7 / ( r_ij^7 + 0.12 sigma_ij^7) ) - 2] *
//*                                                                         *
//***************************************************************************
  double energy,r1,r2,r6,r7,r12,d2;
  double A_term,B_term,AB_term,frac1,frac2,frac3,frac4;
  double sigma3,sigma7,mod;
  sigma3 = sigma*sigma*sigma;
  sigma7 = sigma3*sigma3*sigma;
  VECTOR rij = ri - rj;
  r2 = rij.length2();
  r1 = sqrt(r2);
  r6 = r2 * r2 * r2;
  r7 = r6 * r1;
  frac1 = (1.0/(r7 + 0.12*sigma7));
  B_term = 1.12*sigma7*frac1;

  frac2 = (1.0/(r1 + 0.07*sigma));
  frac3 = 1.07*sigma*frac2;
  frac4 = frac3*frac3;
  A_term = frac4*frac4*frac4*frac3;
  AB_term = A_term*(B_term - 2.0);

  // -dE/dr
  mod = epsilon* 7.0*(frac2*AB_term + r6*frac1*A_term*B_term);
  energy = epsilon * AB_term;

  fi = mod*rij.unit();
  fj = -fj;

}

double Vdw_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                 VECTOR& fi,VECTOR& fj,            /*Outputs*/
                 double D, double r0,double alp){  /*Parameters*/
//********************* double MORSE_STRETCHING *****************************
//*                                                                         *
//*                u = D*{[exp(-alpha*(r_ij-r0))-1]^2 - 1};                 *
//*                                                                         *
//***************************************************************************

  double energy,d;
  VECTOR rij = ri-rj;
  d = rij.length();
  energy = (exp(-alp*(d-r0))-1);

  fi = -energy*D*(-2.0*alp)*(rij/d);
  fj =-fi;

  energy=energy*energy;
  return (energy-1.0)*D;
}


}// namespace libpot
}// liblibra

