#include "RigidBody.h"

namespace librigidbody{

int RigidBody::get_Nf_t(){
  int res = 3;
  if(is_fixed_translation){ res -= 3; }
  return res;
}

int RigidBody::get_Nf_r(){
  int res = 0;
  if(rb_centers_size>=3){ res = 3; }
  else if(rb_centers_size==2){ res = 2; }
  return res;
}

void RigidBody::body_frame_to_lab_frame(const VECTOR& body_frame,VECTOR& lab_frame){
/* This function converts the vector in body frame (time-independent, first)
   into a vector in the lab frame (which center is in point (0,0,0)).
*/
  lab_frame = rb_A_I_to_e_T * body_frame + rb_cm;
}

void RigidBody::lab_frame_to_body_frame(const VECTOR& lab_frame,VECTOR& body_frame){
/* This function converts the vector in lab frame (time-dependent, first)
   into a vector in the body frame (time dependent, second).
*/
  body_frame = - (rb_A_I_to_e * (rb_cm - lab_frame)); //Here is a little trick of VECTOR::operator-
}

VECTOR RigidBody::get_center_in_global_frame(int i){
  return (rb_A_I_to_e_T * rb_centers[i] + rb_cm);
}

VECTOR RigidBody::get_center_in_lab_frame(int i){
  return (rb_A_I_to_e_T * rb_centers[i]);
}

VECTOR RigidBody::get_center_in_body_frame(int i){
  return rb_centers[i];
}

double RigidBody::ekin_rot(){
  return 0.5*(rb_A*rb_l_e.x*rb_l_e.x + rb_B*rb_l_e.y*rb_l_e.y + rb_C*rb_l_e.z*rb_l_e.z);
}

double RigidBody::ekin_tr(){
  return 0.5*rb_iM*rb_p*rb_p;
}


}//namespace librigidbody

