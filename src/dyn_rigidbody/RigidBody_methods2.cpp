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
  \file RigidBody_methods2.cpp
  \brief The file implements RigidBody class interface - modifiers, internal variables transformations/propagation

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{


int RigidBody::scale_angular_(double sc){
/**
  \brief Scale angular variables uniformly by a scalar
  
  This will scale: angular momentum (body frame), angular velocity (body frame), quaternion momentum

  \param[in] sc Scaling parameter
*/
  rb_l_e *= sc;
  rb_w_e *= sc;
  rb_p_r *= sc;
  return 0;
}

int RigidBody::scale_angular_(double scx,double scy,double scz){
/**
  \brief Scale angular variables non-uniformly by scalar components
  
  This will first scale: angular momentum (body frame) and angular velocity (body frame)
  In this case all is more complicated rb_p_r = 2*S(q)*rb_l_e, so rb_p_r is calculated
  directly using scaled rb_l_e variable

  \param[in] scx Scaling parameter for x component of angular velocity and angular momentum
  \param[in] scy Scaling parameter for y component of angular velocity and angular momentum
  \param[in] scz Scaling parameter for z component of angular velocity and angular momentum
*/

  rb_l_e.x *=scx; rb_l_e.y *= scy; rb_l_e.z *= scz;
  rb_w_e.x *=scx; rb_w_e.y *= scy; rb_w_e.z *= scz;

  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::scale_angular_(const MATRIX3x3& sc){
/**
  \brief Scale angular variables non-uniformly by vectors

  We will first transorm angular momentum and velocity (matrix-vector multiplication)
  In this case all is more complicated rb_p_r = 2*S(q)*rb_l_e, so rb_p_r is calculated
  directly using scaled rb_l_e variable

  \param[in] sc Is the transformation (generalized scaling, but actually more than just scaling) matrix.
*/

  rb_l_e = sc*rb_l_e;
  rb_w_e = sc*rb_w_e;

  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0, rb_l_e.x, rb_l_e.y, rb_l_e.z);
  rb_p_r = 2.0*S*tmp;  is_rb_p_r = 1;

  return 0;
}

int RigidBody::scale_linear_(double sc){
/**
  \brief Scale linear variables uniformly by a scalar
  
  This will scale: linear momentum and linear velocity 

  \param[in] sc Scaling parameter
*/

  rb_p *= sc;
  rb_v *= sc;
  return 0;
}

int RigidBody::scale_linear_(double scx,double scy,double scz){
/**
  \brief Scale linear variables non-uniformly by scalar components
  
  This will first scale: linear momentum and linear velocity

  \param[in] scx Scaling parameter for x component of linear velocity and linear momentum
  \param[in] scy Scaling parameter for y component of linear velocity and linear momentum
  \param[in] scz Scaling parameter for z component of linear velocity and linear momentum
*/

  rb_p.x *= scx; rb_p.y *= scy; rb_p.z *= scz;
  rb_v.x *= scx; rb_v.y *= scy; rb_v.z *= scz;
  return 0;
}

int RigidBody::scale_linear_(const MATRIX3x3& sc){
/**
  \brief Scale angular variables non-uniformly by vectors

  We will transorm linear momentum and linear velocity (matrix-vector multiplication)

  \param[in] sc Is the transformation (generalized scaling, but actually more than just scaling) matrix.
*/

  rb_p = sc*rb_p;
  rb_v = sc*rb_v;
  return 0;
}

int RigidBody::scale_position(double sc){
/**
  \brief Scale center of mass coordinates uniformly
  
  \param[in] sc Scaling parameter
*/

  rb_cm *= sc;
  return 0;
}

int RigidBody::scale_position(double scx,double scy,double scz){
/**
  \brief Scale COM non-uniformly by scalar components
  
  \param[in] scx Scaling parameter for x component of COM
  \param[in] scy Scaling parameter for y component of COM
  \param[in] scz Scaling parameter for z component of COM
*/

  rb_cm.x *= scx; rb_cm.y *= scy; rb_cm.z *= scz;
  return 0;
}

int RigidBody::scale_position(const MATRIX3x3& sc){
/**
  \brief Scale COM coordinate by vectors

  We will transorm COM coordinates by a given matrix

  \param[in] sc Is the transformation (generalized scaling, but actually more than just scaling) matrix.
*/

  rb_cm = sc*rb_cm;
  return 0;
}

int RigidBody::shift_angular_momentum(const VECTOR& shft){
/**
  \brief Transform the angular momentum by a vector

  l --> l + shft
  This does not change angular velocity

  \param[in] shft The amount of change
*/

  rb_l_e += shft;
  return 0;
}

int RigidBody::shift_angular_velocity(const VECTOR& shft){
/**
  \brief Transform the angular velocity by a vector

  w --> w + shft
  This does not change angular momentum

  \param[in] shft The amount of change
*/

  rb_w_e += shft;
  return 0;
}

int RigidBody::shift_linear_momentum(const VECTOR& shft){
/**
  \brief Transform the linear momentum by a vector

  p --> p + shft
  This does not change linear velocity

  \param[in] shft The amount of change
*/

  rb_p += shft;
  return 0;
}

int RigidBody::shift_linear_velocity(const VECTOR& shft){
/**
  \brief Transform the linear velocity by a vector

  v --> v + shft
  This does not change linear momentum

  \param[in] shft The amount of change
*/

  rb_v += shft;
  return 0;
}

int RigidBody::shift_position(const VECTOR& shft){
/**
  \brief Transform the COM coordinate by a vector

  r --> r + shft

  \param[in] shft The amount of change
*/

  rb_cm += shft;
  return 0;
}

int RigidBody::fix_translation(){ 
/**
  \brief Freeze linear DOF, so the position does not evolve
*/

  is_fixed_translation = 1; return 0;
}

int RigidBody::fix_rotation(){ 
/**
  \brief Freeze rotational DOF, so the orientation does not evolve
*/

  is_fixed_rotation = 1; return 0;
}

int RigidBody::unfix_translation(){ 
/**
  \brief Un-Freeze linear DOF, so the position can evolve
*/

  is_fixed_translation = 0; return 0;
}

int RigidBody::unfix_rotation(){ 
/**
  \brief Un-Freeze rotational DOF, so the orientation can evolve
*/

  is_fixed_rotation = 0; return 0; 
}


int RigidBody::apply_torque(double dt){
/**
  \brief Apply torque (the one already set in the RB object), to evolve angular momentum

  Angular momentum is in body frame. Quaternion momentum is also evolved. The angular velocity is not evolved.

  l --> l + dt*torque

  \param[in] dt the time of evolution (timestep or fraction of it)
*/

  rb_l_e += dt*rb_torque_e;
  MATRIX S(4,4); init_S_matrix(S);
  QUATERNION tmp(0.0,rb_torque_e.x, rb_torque_e.y, rb_torque_e.z);
  rb_p_r += 2.0*dt*S*tmp;
  return 0;
}

int RigidBody::apply_force(double dt){
/**
  \brief Apply force (the one already set in the RB object), to evolve linear momentum

  The linear velocity is not evolved

  p --> p + dt*f

  \param[in] dt the time of evolution (timestep or fraction of it)
*/

  rb_p += dt*rb_force;
  return 0;
}

int RigidBody::apply_force(MATRIX3x3& x){
/**
  \brief Apply force (the one already set in the RB object), to evolve linear momentum

  The linear velocity is not evolved

  p --> p + x*f

  \param[in] x The scaling factor that is applied to force before it is applied to evolve linear momentum
*/

  rb_p += x*rb_force;
  return 0;
}


}// namespace librigidbody
}// liblibra
