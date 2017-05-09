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

