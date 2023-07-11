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
  \file RigidBody_methods5_4.cpp
  \brief The file implements the KLN propagation method based on 4 in-plane  rotations. 

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{

void RigidBody::propagate_kln(double t){
/**
  \brief The KLN propagation scheme for RB dynamics
  \param[in] t The duration of the propagation

  This function provides an approximate (but not symplectic) solution
  to the free rigid-body problem in terms of 4 consequtive rotations as
  described by:
  Kamberaj, H.; Low, R. J; Neal, M. P. "Time reversible and symplectic
  integrators for molecular dynamics simulations of rigid molecules"
  J. Chem. Phys. 2005, 122, 224114-1 - 224114-30

  Notes on sign in rotation angle:
  direct derivation of action of operator
  iL_i = phi_i *(S_i * d/dS_k - S_k * d/dS_j ) gives:

             S_i = S_i
  exp(iL_i)  S_j = cos(phi_i) * S_j - sin(phi_i) * S_k
             S_k = sin(phi_i) * S_j - cos(phi_i) * S_k

*/
  MATRIX3x3 R;
  VECTOR nx(1.0,0.0,0.0);
  VECTOR ny(0.0,1.0,0.0);

  //--------- 0.5 phiy --------------
//  R.Rotation(0.5*t*rb_l_e.y*(rb_C-rb_B)*ny);
  R.Ry(0.5*t*rb_l_e.y*(rb_C-rb_B));
  rb_l_e      = R * rb_l_e;

  //--------- 0.5 phix --------------
//  R.Rotation(0.5*t*rb_l_e.x*(rb_C-rb_A)*nx);
  R.Rx(0.5*t*rb_l_e.x*(rb_C-rb_A));
  rb_l_e      = R * rb_l_e;

  //---- Quaternion propagation -------
  MATRIX skew1(4,4),mexp(4,4);
  mexp.Init_Unit_Matrix(1.0);

//  set_angular_momentum(rb_l_e);
  rb_w_e.x = rb_l_e.x * rb_A;
  rb_w_e.y = rb_l_e.y * rb_B;
  rb_w_e.z = rb_l_e.z * rb_C;

  double omega = 0.5*t*rb_w_e.length();
  skew1.skew1(rb_w_e);

  // Euler-Rodriguez formula
  mexp = cos(omega)*mexp + 0.5*t*sin_(omega)*skew1;
  rb_L = mexp * rb_L;

  //--------- 0.5 phix --------------
//  R.Rotation(0.5*t*rb_l_e.x*(rb_C-rb_A)*nx);
  R.Rx(0.5*t*rb_l_e.x*(rb_C-rb_A));
  rb_l_e      = R * rb_l_e;

  //--------- 0.5 phiy --------------
//  R.Rotation(0.5*t*rb_l_e.y*(rb_C-rb_B)*ny);
  R.Ry(0.5*t*rb_l_e.y*(rb_C-rb_B));
  rb_l_e      = R * rb_l_e;


  // Update dependent variables
  set_angular_momentum(rb_l_e);
  set_orientation(rb_L);

}

}// namespace librigidbody
}// liblibra
