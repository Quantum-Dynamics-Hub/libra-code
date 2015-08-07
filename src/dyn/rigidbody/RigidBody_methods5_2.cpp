#include "RigidBody.h"

namespace libdyn{
namespace librigidbody{

void RigidBody::propagate_dlml(double t,double& Ps){
/*******************************************************************
 This function provides an approximate (symplectic) solution to the
 free rigid-body problem in terms of 5 consequtive rotations as
 described by:
 Dullweber, A.; Leimkuhler, B.; McLachlan, R.
 "Symplectic splitting methods for rigid body molecular dynamics"
 J. Chem. Phys. 1997, 107, 5840-5851
*******************************************************************/
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
  double Ps = 0.0;
  propagate_dlml(t, Ps);
  return Ps;
}

}// namespace librigidbody
}// namespace libdyn



