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
  \file RigidBody_methods3.cpp
  \brief The file implements RigidBody class interface - getters, properties, external variable transformations

*/

#include "RigidBody.h"

/// liblibra namespace
namespace liblibra{


/// librigidbody namespace
namespace librigidbody{

int RigidBody::get_Nf_t(){
/** Return the number of translational DOF
*/

  int res = 3;
  if(is_fixed_translation){ res -= 3; }
  return res;
}

int RigidBody::get_Nf_r(){
/** Return the number of rotational DOF
*/

  int res = 0;
  if(rb_centers_size>=3){ res = 3; }
  else if(rb_centers_size==2){ res = 2; }
  return res;
}

void RigidBody::body_frame_to_lab_frame(const VECTOR& body_frame,VECTOR& lab_frame){
/** 
  \brief Transformation of a vector from body to lab frame

  This function converts the vector in body frame (time-independent, attached to the principal axes of RB)
  into a vector in the lab frame (which center is in point (0,0,0)).

  \param[in] body_frame The vector in the body frame (the one that moves with the RB)
  \param[out] lab_frame The vector in the lab frame (global, fixed coordinate system)
*/
  lab_frame = rb_A_I_to_e_T * body_frame + rb_cm;
}

void RigidBody::lab_frame_to_body_frame(const VECTOR& lab_frame,VECTOR& body_frame){
/** 
  \brief Transformation of a vector from lab to body frame

  This function converts the vector from lab frame (which center is in point (0,0,0)).
  into a vector in the body frame (time-independent, attached to the principal axes of RB)

  \param[in] lab_frame The vector in the lab frame (global, fixed coordinate system)
  \param[out] body_frame The vector in the body frame (the one that moves with the RB)

*/
  body_frame = - (rb_A_I_to_e * (rb_cm - lab_frame)); //Here is a little trick of VECTOR::operator-
}

VECTOR RigidBody::get_center_in_global_frame(int i){
/** 
  \brief Return the coordinates (Cartesian) of a material point center - in the global lab frame (center at point 0,0,0)

  \param[in] i The index of the material point (center) included in RB. For this center, we return the global lab frame coordinate

*/

  return (rb_A_I_to_e_T * rb_centers[i] + rb_cm);
}

VECTOR RigidBody::get_center_in_lab_frame(int i){
/** 
  \brief Return the coordinates (Cartesian) of a material point center - in the moving lab frame (center is at the COM of RB)

  \param[in] i The index of the material point (center) included in RB. For this center, we return the moving lab frame coordinate

*/

  return (rb_A_I_to_e_T * rb_centers[i]);
}

VECTOR RigidBody::get_center_in_body_frame(int i){
/** 
  \brief Return the coordinates (Cartesian) of a material point center - in the body frame

  This coordinate is fixed for the point in the RB. It doesn't change when RB moves and rotates.

  \param[in] i The index of the material point (center) included in RB. For this center, we return the body frame coordinate

*/

  return rb_centers[i];
}

double RigidBody::ekin_rot(){
/** 
  \brief Return the rotational kinetic energy of RB

*/

  return 0.5*(rb_A*rb_l_e.x*rb_l_e.x + rb_B*rb_l_e.y*rb_l_e.y + rb_C*rb_l_e.z*rb_l_e.z);
}

double RigidBody::ekin_tr(){
/** 
  \brief Return the translational kinetic energy of RB

*/

  return 0.5*rb_iM*rb_p*rb_p;
}


}// namespace librigidbody
}// liblibra
