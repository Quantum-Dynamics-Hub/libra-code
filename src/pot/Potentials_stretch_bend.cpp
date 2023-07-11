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

#include "Potentials_stretch_bend.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{

double Stretch_Bend_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3,         /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3,         /* Outputs*/
                             double k_ijk,double k_kji, double theta_0,
                             double r_ij0,double r_kj0                 /* Parameters*/
                            ){
//******************** double MMFF94 type STRETCH-BEND TERM **********************
//*                                                                              *
//*            E = 2.51210*(theta_ijk-theta_0)*[k_ijk *drij + k_kji*dr_kj];      *
//*                                                                              *
//*  theta_0 - in radians                                                        *
//*                                                                              *
//********************************************************************************
  VECTOR r12 = r1 - r2;
  VECTOR r32 = r3 - r2;
  VECTOR f12,f32;
  double d12,d32,cos_theta,sin_theta,theta,diff,diff1,diff3;
  double energy;
  d12 = r12.length();
  d32 = r32.length();

  // Calculate angle
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta >  1.0) { cos_theta = 1.0;}
  else if (cos_theta < -1.0) { cos_theta = -1.0;}
  sin_theta = sqrt(1.0 - cos_theta * cos_theta);
  theta = acos(cos_theta);

  // Energy
  diff = rad_to_deg*(theta-theta_0);
  diff1 = 2.51210  * k_ijk * diff;
  diff3 = 2.51210  * k_kji * diff;
  energy = diff1*(d12 - r_ij0) + diff3*(d32 - r_kj0);

  // These are:
  // - dE/dcos(theta) *(dcos(theta)/drx) x = i,j,k
  // and simply -dE/dr terms
  f1 = - diff1 * r12.unit();
  diff1 = (-1.0/sin_theta)* 2.51210  * k_ijk * (d12 - r_ij0);
  f12 = (diff1/d12)*(r12.unit()*cos_theta - r32.unit());
  f1 += f12;

  f3 = - diff3 * r32.unit();
  diff3 = (-1.0/sin_theta)* 2.51210  * k_kji * (d32 - r_kj0);
  f32 = (diff/d32)*(r32.unit()*cos_theta - r12.unit());
  f3 += f32;

  f2 = -f12 - f32;

  return energy;
}


}//namespace libpot
}// liblibra

