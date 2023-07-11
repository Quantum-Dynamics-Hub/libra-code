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

#include "Switching_functions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{

void SWITCH(VECTOR& r1,VECTOR& r2,     /*Input*/
            double R_on,double R_off,  /*Paramaters*/
            double& SW,VECTOR& dSW){   /*Output*/
/*****************************************************************************
  E = ((R_off-R)/(R_off - R_on))^3 *
      ( 1 + 3*(R-R_on)/(R_off - R_on) + 6*((R-R_on)/(R_off - R_on))^2 )
*****************************************************************************/
  VECTOR r12 = r1 - r2;
  double dist1 = r12.length();
  if(dist1<=R_on){ SW = 1.0; dSW = 0.0;   }
  else if((dist1>R_on)&&(dist1<R_off)){
    double dR12 = (R_off - R_on);
    double x = ((R_off-dist1)/dR12);
    double y = ((dist1-R_on)/dR12);
    double x2 = x*x;
    double Y = (1.0 + 3.0*y + 6.0*y*y);

    SW = x2*x*Y;
    dSW = 3.0*(x2/dR12)*( x*(1.0 + 4.0*y) - Y )*(r12/dist1);
  }
  else if(dist1>=R_off){ SW = 0.0;   dSW = 0.0;    }

}

void DOUBLE_SWITCH(double x,double a,double eps,double& SW,double& dSW){
  double ksi1,ksi2;
  if((x>=0.0)&&(x<=eps)){
    ksi1 = (x/eps);
    ksi2 = (1.0 - ksi1);
    SW = ksi1*ksi1*ksi1*(1.0 + 3.0*ksi2 + 6.0*ksi2*ksi2);
    dSW = (30.0/eps)*ksi1*ksi1*ksi2*ksi2;
  }
  else if((x>eps)&&(x<(a-eps))){  SW = 1.0;  dSW = 0.0;   }
  else if((x>=(a-eps))&&(x<=a)){
    ksi1 = ((a-x)/eps);
    ksi2 = (1.0 - ksi1);
    SW = ksi1*ksi1*ksi1*(1.0 + 3.0*ksi2 + 6.0*ksi2*ksi2);
    dSW = -(30.0/eps)*ksi1*ksi1*ksi2*ksi2;
  }
  else{ SW = 0.0;  dSW = 0.0;   }

}



// Here go just wrapper functions for export into Python
boost::python::list SWITCH(VECTOR r1,VECTOR r2, double R_on,double R_off){

  boost::python::list res;
  double SW = 0.0;
  VECTOR dSW;

  SWITCH(r1,r2,R_on,R_off,SW,dSW);

  res.append(SW);
  res.append(dSW);
 
  return res;

}


boost::python::list DOUBLE_SWITCH(double x,double a,double eps){

  boost::python::list res;
  double SW = 0.0;
  double dSW = 0.0;

  DOUBLE_SWITCH(x, a, eps,SW,dSW);

  res.append(SW);
  res.append(dSW);
 
  return res;

}




}// namespace libpot
}// liblibra

