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
  \file RigidBody_methods5_6.cpp
  \brief The file implements the method of Omelyan for RB propagation 
  
*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{

void RigidBody::propagate_omelyan(double t){
/**
  \brief Omelyan integrator
  \param[in] t Integration duration

  This function provides an approximate solution to the torqued
  rigid-body problem in terms of leap-frog-like algorithm  as
  described by:
  Omelyan I. P. "Algorithm for numerical integration of the rigid-body
  equations of motion" Phys. Rev. E. 1998, V. 58, P. 1169-1172
  Here we, however, consider a free rigid-body propagation
  The variables propagate here according to the following scheme:
  l_e (w_e):  t-h/2 -> t+h/2
  A_I_to_e :  t     -> t+h
*/

//  set_angular_momentum(rb_l_e);
// l_e = l_e(t-h/2)
  rb_w_e.x = rb_l_e.x * rb_A;
  rb_w_e.y = rb_l_e.y * rb_B;
  rb_w_e.z = rb_l_e.z * rb_C;


  MATRIX3x3 D,W,P;
  VECTOR omega0(rb_w_e),omega1,omega_e(rb_w_e);
  double om2,mod1,mod2,err;
  double Jx = rb_I_e.xx;
  double Jy = rb_I_e.yy;
  double Jz = rb_I_e.zz;

  // ----- Leap-frog-like iteration procedure to obtain angular velocities on half step--------
  do{
    omega1.x = omega_e.x + rb_A*t*(rb_torque_e.x + 0.5*(Jy - Jz)*(omega_e.y*omega_e.z + omega0.y*omega0.z) );
    omega1.y = omega_e.y + rb_B*t*(rb_torque_e.y + 0.5*(Jz - Jx)*(omega_e.x*omega_e.z + omega0.x*omega0.z) );
    omega1.z = omega_e.z + rb_C*t*(rb_torque_e.z + 0.5*(Jx - Jy)*(omega_e.x*omega_e.y + omega0.x*omega0.y) );
    err = (omega1-omega0).length2();
    omega0 = omega1;
  }while(err>MACHPREC);
  omega_e = omega1;
  rb_w_e = omega_e;

  //---- Propagate of orientational variables -------------
  om2 = omega_e.length2();
  mod1 = (1.0 - 0.25*t*t*om2);
  mod2 = (1.0 + 0.25*t*t*om2);
  D.diag(mod1);
  W.skew(-omega_e);
  P.tensor_product(omega_e,omega_e);
  D = D + t*(W + 0.5*t*P);
  D = D/mod2;
  // rb_A_I_to_e: t -> t+h
  rb_A_I_to_e = D*rb_A_I_to_e;

  // Update dependent variables
//  set_angular_velocity(rb_w_e);
  rb_l_e.x = rb_w_e.x / rb_A;
  rb_l_e.y = rb_w_e.y / rb_B;
  rb_l_e.z = rb_w_e.z / rb_C;

  set_orientation(rb_A_I_to_e);

}

}// namespace librigidbody
}// liblibra
