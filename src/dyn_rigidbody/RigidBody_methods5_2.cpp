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
/**
  \file RigidBody_methods5_2.cpp
  \brief The file implements the DLML propagation method based on 5 rotations 

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{

void RigidBody::propagate_dlml(double t,double& Ps){
/** 
  \brief DLML propagation scheme
  \param[in] t Propagation duration
  \param[out] Ps The variable that is needed for Nose-Poincare thermostat evolvution, when 
                 DLML is coupled to it   

  This function provides an approximate (symplectic) solution to the
  free rigid-body problem in terms of 5 consequtive rotations as
  described by:
  Dullweber, A.; Leimkuhler, B.; McLachlan, R.
  "Symplectic splitting methods for rigid body molecular dynamics"
  J. Chem. Phys. 1997, 107, 5840-5851
*/
  MATRIX3x3 R;
  VECTOR nx(1.0,0.0,0.0);
  VECTOR ny(0.0,1.0,0.0);
  VECTOR nz(0.0,0.0,1.0);
  Ps = 0.0;
  double phi;

  // This version is based on l_e!!!

  //--------- 0.5 phix --------------
  phi = -0.5*t*rb_A*rb_l_e.x;
  Ps -= phi*(0.5*rb_l_e.x);
  R.Rx(phi);
  rb_l_e      = R * rb_l_e;
  rb_A_I_to_e = R * rb_A_I_to_e;

  //--------- 0.5 phiy --------------
  phi = -0.5*t*rb_B*rb_l_e.y;
  Ps -= phi*(0.5*rb_l_e.y);
  R.Ry(phi);
  rb_l_e      = R * rb_l_e;
  rb_A_I_to_e = R * rb_A_I_to_e;

  //--------- 1.0 phiz --------------
  phi = -t*rb_C*rb_l_e.z;
  Ps -= phi*(0.5*rb_l_e.z);
  R.Rz(phi);
  rb_l_e      = R * rb_l_e;
  rb_A_I_to_e = R * rb_A_I_to_e;

  //--------- 0.5 phiy --------------
  phi = -0.5*t*rb_B*rb_l_e.y;
  Ps -= phi*(0.5*rb_l_e.y);
  R.Ry(phi);
  rb_l_e      = R * rb_l_e;
  rb_A_I_to_e = R * rb_A_I_to_e;

  //--------- 0.5 phix --------------
  phi = -0.5*t*rb_A*rb_l_e.x;
  Ps -= phi*(0.5*rb_l_e.x);
  R.Rx(phi);
  rb_l_e      = R * rb_l_e;
  rb_A_I_to_e = R * rb_A_I_to_e;


  // Update dependent variables
//  set_angular_velocity(rb_w_e);
  set_angular_momentum(rb_l_e);
  set_orientation(rb_A_I_to_e);

}

double RigidBody::propagate_dlml(double t){
/** 
  \brief DLML propagation scheme
  \param[in] t Propagation duration

  This function provides an approximate (symplectic) solution to the
  free rigid-body problem in terms of 5 consequtive rotations as
  described by:
  Dullweber, A.; Leimkuhler, B.; McLachlan, R.
  "Symplectic splitting methods for rigid body molecular dynamics"
  J. Chem. Phys. 1997, 107, 5840-5851
*/

  double Ps = 0.0;
  propagate_dlml(t, Ps);
  return Ps;
}

}// namespace librigidbody

}// liblibra

