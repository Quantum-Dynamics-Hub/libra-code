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

#include "Potentials_bonds.h"


/// liblibra namespace
namespace liblibra{


namespace libpot{


double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj,  /*Outputs*/
                      double K, double r0){  /*Parameters*/
//****************** double HARMONIC_STRETCHING ******************************
//*                                                                          *
//*       E =           K*(r_ij-r0)^2;                                       *
//*                                                                          *
//****************************************************************************

  double energy,d;
  VECTOR rij = ri - rj;
  d = rij.length();
  energy = (d-r0);

  fi = -2.0*K*energy*(rij/d);
  fj = -fi;

  return K*energy*energy;
}

double Bond_Harmonic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                     VECTOR& fi,VECTOR& fj, MATRIX& Hess,  /*Outputs*/
                      double K, double r0, int opt){  /*Parameters*/
//****************** double HARMONIC_STRETCHING *********************************
//*                                                                             *
//*       E =           K*(|r_ij|-r0)^2;                                        *
//*                                                                             *
//*   u_alp = (r_ij)_alp/|r_ij|                                                 *
//*                                                                             *
//*   du_alp/dri_bet = -(rij)_alp * (rij)_bet/|r_ij|^3  + 1/|r_ij| * delta_alp_bet  *  
//*                                                                                 *
//*  =  (delta_alp_bet - u_alp * u_bet )/|r_ij|                                     *
//*                                                                                 *
//*   du_alp/drj_bet =  (rij)_alp * (rij)_bet/|r_ij|^3  - 1/|r_ij| * delta_alp_bet  *  
//*                                                                                 *
//*  = -(delta_alp_bet - u_alp * u_bet )/|r_ij|                                     *
//*                                                                                 *
//*                                                                             *
//*  dE_ri_alp = 2*K*(|r_ij| - r0) * u_alp                                      *
//*  dE_rj_alp =-2*K*(|r_ij| - r0) * u_alp                                      *
//*                                                                             *
//*  dE_ri_alp_dri_bet = 2*K*[ u_alp * u_bet + (|r_ij| - r0) * du_alp/dri_bet ] *
//*  dE_ri_alp_drj_bet = 2*K*[-u_alp * u_bet + (|r_ij| - r0) * du_alp/drj_bet ] *
//*                                                                             *
//*******************************************************************************

  double energy,d, dr;
  VECTOR rij = ri - rj;
  d = rij.length();
  dr = (d-r0);

  // Energy
  energy = 0.0;
  
  if(opt>=0){

    energy = K*dr*dr;


    // Forces

    if(opt>=1){

      VECTOR u = rij/d; 
 
      fi = -2.0*K*dr*u;   // -dE/dr_i
      fj = -fi;           // -dE/dr_j


        //          d/dxi   d/dyi    d/dzi      d/dxj    d/dyj    d/dzj
        // d/dxi
        // d/dyi 
        // d/dzi 
        // d/dxj                       Hessian
        // d/dyj
        // d/dzj
        
      
      if(opt>=2){

        
        Hess.M[0] = u.x * u.x + dr * (1.0 - u.x * u.x)/d;   // d^2 E / dxi * dxi
        Hess.M[1] = u.x * u.y - dr * u.x * u.y/d;           // d^2 E / dxi * dyi
        Hess.M[2] = u.x * u.z - dr * u.x * u.z/d;           // d^2 E / dxi * dzi
        Hess.M[3] = -Hess.M[0];                             // d^2 E / dxi * dxj
        Hess.M[4] = -Hess.M[1];                             // d^2 E / dxi * dyj
        Hess.M[5] = -Hess.M[2];                             // d^2 E / dxi * dzj
        
        Hess.M[6] = u.y * u.x - dr * u.y * u.x/d;           // d^2 E / dyi * dxi
        Hess.M[7] = u.y * u.y - dr * (1.0 - u.y * u.y)/d;   // d^2 E / dyi * dyi
        Hess.M[8] = u.y * u.z - dr * u.y * u.z/d;           // d^2 E / dyi * dzi
        Hess.M[9] = -Hess.M[6];                             // d^2 E / dyi * dxj
        Hess.M[10]= -Hess.M[7];                             // d^2 E / dyi * dyj
        Hess.M[11]= -Hess.M[8];                             // d^2 E / dyi * dzj
        
        Hess.M[12]= u.z * u.x - dr * u.z * u.x/d;           // d^2 E / dzi * dxi
        Hess.M[13]= u.z * u.y - dr * u.z * u.y/d;           // d^2 E / dzi * dyi
        Hess.M[14]= u.z * u.z + dr * (1.0 - u.z * u.z)/d;   // d^2 E / dzi * dzi
        Hess.M[15]= -Hess.M[12];                            // d^2 E / dyi * dxj
        Hess.M[16]= -Hess.M[13];                            // d^2 E / dyi * dyj
        Hess.M[17]= -Hess.M[14];                            // d^2 E / dyi * dzj

        for(int i=18;i<36;i++){  Hess.M[i] = -Hess.M[i-18]; } // by symmetry

 
        Hess *= 2.0*K;  

      }// opt >= 2
    }// opt >= 1
  }// opt >= 0

  return energy;
}


double Bond_Quartic(VECTOR& ri,VECTOR& rj,  /*Inputs*/
                    VECTOR& fi,VECTOR& fj,  /*Outputs*/
                    double K, double r0){   /*Parameters*/
//***************** double QUARTIC BOND STRETCH *****************************
//*                                                                         *
//*                u = K*(r-r0)^2*[1+cs*(r-r0)+7/12*cs^2*(r-r0)^2]          *
//*                                                                         *
//* Used in: MMFF94                                                         *
//***************************************************************************
  double cs = -2.0;
  double cs2 = (7.0/12.0)*cs*cs;
  double d,d1,d2;
  VECTOR rij = ri-rj;
  d = rij.length();
  d1 = d - r0;
  d2 = d1*d1;

  fi = -K*d*(2.0 + 3.0*cs*d + 4.0*cs2*d2)*(rij/d);
  fj = -fi;

  return K*d2*(1.0 + cs*d + cs2*d2);
}

double Bond_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                  VECTOR& fi,VECTOR& fj,            /*Outputs*/
                  double D, double r0,double alp){  /*Parameters*/
//********************* double MORSE_STRETCHING *****************************
//*                                                                         *
//*                u = D*{[exp(-alpha*(r_ij-r0))-1]^2 - 1};                 *
//*                                                                         *
//*                u = D*{x^2 - 2.0*x},where x = exp(-alpha*(r_ij-r0))      *
//*                                                                         *
//***************************************************************************

  double energy,d;
  VECTOR rij = ri-rj;
  d = rij.length();
/*
  energy = (exp(-alp*(d-r0))-1);

  fi = -energy*D*(-2.0*alp)*(rij/d);
  fj =-fi;

  energy=energy*energy;
  return (energy-1.0)*D;
*/
  double x = std::exp(-alp*(d-r0));

  fi = 2.0*D*(x - 1.0)*(alp/d)*rij;
  fj =-fi;
  
  return D*(x - 2.0)*x;
}





boost::python::list Bond_Harmonic(VECTOR ri,VECTOR rj, double K, double r0){

  boost::python::list res;
  double en = 0.0;
  VECTOR fi, fj;

  en = Bond_Harmonic(ri,rj,fi,fj,K,r0);

  res.append(en);
  res.append(fi);
  res.append(fj);
 
  return res;

}


boost::python::list Bond_Harmonic(VECTOR ri,VECTOR rj, double K, double r0, int opt){
/**
  opt = 0 - energy 
*/

  boost::python::list res;
  double en = 0.0;
  VECTOR fi, fj;
  MATRIX Hess(6,6);
  
  en = Bond_Harmonic(ri,rj,fi,fj,Hess,K,r0,opt);

  res.append(en);
  res.append(fi);
  res.append(fj);
  res.append(Hess);
 
  return res;

}


boost::python::list Bond_Quartic(VECTOR ri,VECTOR rj, double K, double r0){

  boost::python::list res;
  double en = 0.0;
  VECTOR fi, fj;

  en = Bond_Quartic(ri,rj,fi,fj,K,r0);

  res.append(en);
  res.append(fi);
  res.append(fj);
 
  return res;

}

boost::python::list Bond_Morse(VECTOR ri,VECTOR rj, double D, double r0,double alp){

  boost::python::list res;
  double en = 0.0;
  VECTOR fi, fj;

  en = Bond_Morse(ri,rj,fi,fj,D,r0,alp);

  res.append(en);
  res.append(fi);
  res.append(fj);
 
  return res;

}




}// namespace libpot
}// liblibra

