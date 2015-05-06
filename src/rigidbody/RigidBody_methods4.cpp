#include "RigidBody.h"

namespace librigidbody{

void RigidBody::Rotate_I_x(double phi){
/************************************************
 Rotation around center of mass along the x
 axis in external coordinate system (lab frame)
************************************************/
  MATRIX3x3 R; R.Rx(phi);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate_I_y(double phi){
/************************************************
 Rotation around center of mass along the y
 axis in external coordinate system (lab frame)
************************************************/
  MATRIX3x3 R; R.Rx(phi);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate_I_z(double phi){
/************************************************
 Rotation around center of mass along the z
 axis in external coordinate system (lab frame)
************************************************/
  MATRIX3x3 R; R.Rx(phi);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate_e_x(double phi){
/************************************************
 Rotation around center of mass along the x
 axis in body-fixed coordinate system (body frame)
************************************************/
  VECTOR ux,uy,uz;
  rb_A_I_to_e.get_vectors(ux,uy,uz);
  MATRIX3x3 R; R.Rotation(phi*ux);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate_e_y(double phi){
/************************************************
 Rotation around center of mass along the y
 axis in body-fixed coordinate system (body frame)
************************************************/
  VECTOR ux,uy,uz;
  rb_A_I_to_e.get_vectors(ux,uy,uz);
  MATRIX3x3 R; R.Rotation(phi*ux);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate_e_z(double phi){
/************************************************
 Rotation around center of mass along the y
 axis in body-fixed coordinate system (body frame)
************************************************/
  VECTOR ux,uy,uz;
  rb_A_I_to_e.get_vectors(ux,uy,uz);
  MATRIX3x3 R; R.Rotation(phi*ux);
  rb_A_I_to_e = R * rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate(MATRIX3x3& R){
/*************************************************
 Warning: Matrix R should be a valid rotation matrix
*************************************************/
  rb_A_I_to_e = R*rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate(QUATERNION& quat){
/************************************************
 Rotation of rigid body by quaternion quat
*************************************************/
  MATRIX3x3 R;
  QUATERNION_TO_MATRIX(quat,R);
  rb_A_I_to_e = R*rb_A_I_to_e;
  set_orientation(rb_A_I_to_e);
}

void RigidBody::Rotate(double phi,VECTOR& dir){
/*************************************************
 Rotation by angle phi around direction dir
**************************************************/
  double cs = cos(0.5*phi);
  double si = sin(0.5*phi);
  VECTOR u = dir.unit();
  QUATERNION quat(cs,si*u.x,si*u.y,si*u.z);
  Rotate(quat);
}


}//namespace librigidbody

