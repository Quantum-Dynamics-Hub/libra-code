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

#include "RigidBody.h"

namespace libdyn{
namespace librigidbody{


int RigidBody::scale_angular_(double sc){
  rb_l_e *= sc;
  rb_w_e *= sc;
  rb_p_r *= sc;
  return 0;
}

int RigidBody::scale_angular_(double scx,double scy,double scz){
  rb_l_e.x *=scx; rb_l_e.y *= scy; rb_l_e.z *= scz;
  rb_w_e.x *=scx; rb_w_e.y *= scy; rb_w_e.z *= scz;
  // In this case all is more complicated p_r = 2*S(q)*l_e, p_r is calculated
  // directly using scaled l_e variables:
  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::scale_angular_(const MATRIX3x3& sc){
  rb_l_e = sc*rb_l_e;
  rb_w_e = sc*rb_w_e;
  // In this case all is more complicated p_r = 2*S(q)*l_e, p_r is calculated
  // directly using scaled l_e variables:
  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::scale_linear_(double sc){
  rb_p *= sc;
  rb_v *= sc;
  return 0;
}

int RigidBody::scale_linear_(double scx,double scy,double scz){
  rb_p.x *= scx; rb_p.y *= scy; rb_p.z *= scz;
  rb_v.x *= scx; rb_v.y *= scy; rb_v.z *= scz;
  return 0;
}

int RigidBody::scale_linear_(const MATRIX3x3& sc){
  rb_p = sc*rb_p;
  rb_v = sc*rb_v;
  return 0;
}

int RigidBody::scale_position(double sc){
  rb_cm *= sc;
  return 0;
}

int RigidBody::scale_position(double scx,double scy,double scz){
  rb_cm.x *= scx; rb_cm.y *= scy; rb_cm.z *= scz;
  return 0;
}

int RigidBody::scale_position(const MATRIX3x3& sc){
  rb_cm = sc*rb_cm;
  return 0;
}

int RigidBody::shift_angular_momentum(const VECTOR& shft){
  rb_l_e += shft;
  return 0;
}

int RigidBody::shift_angular_velocity(const VECTOR& shft){
  rb_w_e += shft;
  return 0;
}

int RigidBody::shift_linear_momentum(const VECTOR& shft){
  rb_p += shft;
  return 0;
}

int RigidBody::shift_linear_velocity(const VECTOR& shft){
  rb_v += shft;
  return 0;
}

int RigidBody::shift_position(const VECTOR& shft){
  rb_cm += shft;
  return 0;
}

int RigidBody::fix_translation(){ is_fixed_translation = 1; return 0;}
int RigidBody::fix_rotation(){ is_fixed_rotation = 1; return 0;}
int RigidBody::unfix_translation(){ is_fixed_translation = 0; return 0;}
int RigidBody::unfix_rotation(){ is_fixed_rotation = 0; return 0; }


int RigidBody::apply_torque(double dt){
  rb_l_e += dt*rb_torque_e;
  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0,rb_torque_e.x, rb_torque_e.y, rb_torque_e.z);
  rb_p_r += 2.0*dt*S*tmp;
  return 0;
}

int RigidBody::apply_force(double dt){
  rb_p += dt*rb_force;
  return 0;
}

int RigidBody::apply_force(MATRIX3x3& x){
  rb_p += x*rb_force;
  return 0;
}


}// namespace librigidbody
}// namespace libdyn

